"""Inference-only extraction of Sanofi RiboNN. See MODIFICATIONS.md and LICENSE."""
import torch
from torch import nn


class ConvBlock(nn.Module):
    """Conv1D block with or without residual layers"""

    def __init__(
        self,
        in_channels: int,
        out_channels: int,
        kernel_size: int,
        stride: int,
        padding: int,
        eps: float,
        dropout: float,
        residual: bool = False,
        is_initial: bool = False,
        activation=nn.ReLU,
    ):
        super().__init__()
        self.residual = residual
        self.is_initial = is_initial
        self.activation = activation

        if self.is_initial:
            self.conv = nn.Conv1d(
                in_channels=in_channels,
                out_channels=out_channels,
                kernel_size=kernel_size,
                padding=padding,
                stride=stride,
                bias=False,
            )
            if self.residual:
                # (N, L, C) --> (N, L, C), normalize across the last dimension
                self.residual_layernorm = nn.LayerNorm(
                    normalized_shape=out_channels, eps=eps
                )
                self.residual_conv = nn.Sequential(
                    self.activation(),
                    nn.Conv1d(
                        in_channels=out_channels,
                        out_channels=out_channels,
                        kernel_size=3,
                        stride=stride,
                        padding="same",
                    ),
                    nn.Dropout(dropout),
                )
                self.scale = nn.Parameter(torch.zeros(out_channels))
        else:
            # (N, L, C) --> (N, L, C), normalize across the last dimension
            self.layernorm = nn.LayerNorm(normalized_shape=in_channels, eps=eps)
            self.conv = nn.Sequential(
                self.activation(),
                nn.Conv1d(
                    in_channels=in_channels,
                    out_channels=out_channels,
                    kernel_size=kernel_size,
                    stride=stride,
                    padding=padding,
                ),
                nn.Dropout(dropout),
            )
            self.pool = nn.MaxPool1d(kernel_size=2, stride=2)

            if self.residual:
                # (N, L, C) --> (N, L, C), normalize across the last dimension
                self.residual_layernorm = nn.LayerNorm(
                    normalized_shape=out_channels, eps=eps
                )
                self.residual_conv = nn.Sequential(
                    self.activation(),
                    nn.Conv1d(
                        in_channels=out_channels,
                        out_channels=out_channels,
                        kernel_size=3,
                        stride=stride,
                        padding="same",
                    ),
                    nn.Dropout(dropout),
                )
                self.scale = nn.Parameter(torch.zeros(out_channels))

    def forward(self, x):
        """The forward pass

        Args:
            x (Tensor of shape: (N, C, L1)): input tensor

        Returns:
            Tensor of shape: (N, C, L2)
        """
        if self.is_initial:
            # (N, C, L1) --> (N, C, L2)
            x = self.conv(x)

            if self.residual:
                # (N, C, L) --> (N, L, C)
                xt = x.transpose(1, 2)
                # (N, L, C) --> (N, L, C)
                xt = self.residual_layernorm(xt)
                # (N, L, C) --> (N, C, L)
                xt = xt.transpose(1, 2)
                # scale in the channel dimension
                x += self.residual_conv(xt) * self.scale[None, :, None]

            return x

        # (N, C, L) --> (N, L, C)
        x = x.transpose(1, 2)
        # (N, L, C) --> (N, L, C)
        x = self.layernorm(x)
        # (N, L, C) --> (N, C, L)
        x = x.transpose(1, 2)
        # (N, C, L1) --> (N, C, L2)
        x = self.conv(x)

        if self.residual:
            # (N, C, L) --> (N, L, C)
            xt = x.transpose(1, 2)
            # (N, L, C) --> (N, L, C)
            xt = self.residual_layernorm(xt)
            # (N, L, C) --> (N, C, L)
            xt = xt.transpose(1, 2)
            # scale in the channel dimension
            x += self.residual_conv(xt) * self.scale[None, :, None]

        return self.pool(x)


class RiboNN(nn.Module):
    """Frozen upstream convolutional architecture for inference only."""

    def __init__(self, **hparams):
        super().__init__()
        self.activation = nn.ReLU
        for k, v in hparams.items():
            if k == "activation":
                if v.lower() in ("silu", "swish"):
                    v = nn.SiLU
                elif v.lower() == "leakyrelu":
                    v = nn.LeakyReLU
                else:
                    continue
            setattr(self, k, v)

        if self.max_shift != 0:
            raise ValueError("Frozen inference requires max_shift=0")

        # number of channels of initial convolution
        seq_channels = 4 * 3 if self.split_utr5_cds_utr3_channels else 4
        label_channels = (
            self.label_codons
            + self.label_utr5
            + self.label_utr3
            + self.label_splice_sites
            + self.label_up_probs
        )
        input_channels = seq_channels + label_channels

        # shared initial convolution
        self.initial_conv = ConvBlock(
            in_channels=input_channels,
            out_channels=self.filters,
            kernel_size=5,
            stride=self.conv_stride,
            padding=self.conv_padding,
            eps=self.ln_epsilon,
            dropout=self.dropout,
            residual=self.residual,
            is_initial=True,
            activation=self.activation,
        )

        # shared middle convolutions
        self.middle_convs = nn.ModuleList(
            [
                ConvBlock(
                    in_channels=self.filters,
                    out_channels=self.filters,
                    kernel_size=self.kernel_size,
                    stride=self.conv_stride,
                    padding=self.conv_padding,
                    eps=self.ln_epsilon,
                    dropout=self.dropout,
                    residual=self.residual,
                    is_initial=False,
                    activation=self.activation,
                )
                for i in range(self.num_conv_layers)
            ]
        )

        self.head = nn.Sequential(
            self.activation(),
            # nn.AdaptiveMaxPool1d(8),
            nn.Flatten(),
            nn.Dropout(self.dropout),
            nn.Linear(self.filters * self.len_after_conv, self.filters, bias=False),
            self.activation(),
            nn.BatchNorm1d(self.filters),
            nn.Dropout(self.dropout),
            nn.Linear(self.filters, self.num_targets),
        )

    def forward(self, x):
        """The forward pass

        Args:
            x: input tensor of shape: (N, C, L)
            task_idx: task index

        Returns:
            Tensor of shape (N, self.num_targets)
        """
        # augmentation
        # (N, C, L) --> (N, C, L)

        # initial convolutions
        # (N, C1, L1) --> (N, C2, L2)
        x = self.initial_conv(x)

        # middle convolutions
        # (N, C1, L1) --> (N, C2, L2)
        for l in self.middle_convs:
            x = l(x)

        x = self.head(x)

        return x
