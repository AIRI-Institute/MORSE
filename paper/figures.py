"""Rebuild every scientific panel from the frozen individual numerical records.

Figure 1 is a static conceptual vector source. Numerical Figures 2–5 and
Extended Data 1–3 never read a previous PDF. Their layout is reconstructed;
the exact published PDFs remain immutable visual references.
"""
from __future__ import annotations

from collections import defaultdict
import json
from pathlib import Path
import shutil

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
from scipy.stats import rankdata, spearmanr

from .io import DATA, HERE, read_json, rows, runs, sha256, write_json

INK = "#172333"
COLORS = ["#9870B1", "#509893", "#BC8A55", "#4395B5", "#507FB1", "#A64D79", "#8D9355"]
LABELS = {"initial": "Initial", "mfe_per_nt": "MFE", "translation_efficiency_score": "Orthrus TE", "csi_percentile": "CSI",
          "ribodecode_translation_log": "RiboDecode", "ribonn_mean_te": "RiboNN", "best_hv": "Best HV", "maximin": "Maximin", "cai_log": "CAI"}
NAMES = ["Figure1_design_workflow", "Figure2_Mauger_U_rank_heatmap", "Figure3_modular_design", "Figure4_comparative_design",
         "Figure5_predictor_disagreement", "ExtendedDataFigure1_feasibility", "ExtendedDataFigure2_LinearDesign_UTR_control", "ExtendedDataFigure3_RiboDecode_pairs"]


def style():
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7, "axes.labelsize": 7,
                         "axes.titlesize": 8, "legend.fontsize": 6, "xtick.labelsize": 6, "ytick.labelsize": 6,
                         "text.color": INK, "axes.labelcolor": INK, "pdf.fonttype": 42,
                         "axes.spines.top": False, "axes.spines.right": False, "savefig.facecolor": "white"})


def title(ax, letter, text):
    ax.set_title(f"{letter}  {text}", loc="left", fontweight="bold", pad=10)


def save(fig, output, name, fmt):
    if fmt in {"pdf", "both"}:
        fig.savefig(output / (name + ".pdf"), bbox_inches="tight", pad_inches=.12, metadata={"CreationDate": None, "ModDate": None})
    if fmt in {"png", "both"}:
        fig.savefig(output / (name + ".png"), dpi=220, bbox_inches="tight", pad_inches=.12)
    plt.close(fig)


def scope(ax, mutable, scorer_rows, *, name=None):
    parts = ["5′UTR", "AUG", "CDS", "UAA", "3′UTR", "A100"]
    widths = np.asarray([1.2, .55, 3, .55, 1.2, .7])
    x = np.r_[0, np.cumsum(widths)]
    for i, part in enumerate(parts):
        ax.add_patch(Rectangle((x[i], 0), widths[i]-.05, .44, facecolor="#EDE9FE" if i in mutable else "#EDF1F5", edgecolor="#B8C3CD", lw=.5))
        ax.text(x[i]+(widths[i]-.05)/2, .22, part, ha="center", va="center", fontsize=7)
    for i, (label, first, last, color) in enumerate(scorer_rows):
        y = -.42-i*.35
        ax.plot([x[first], x[last]], [y, y], lw=2, color=color)
        ax.text(x[-1]+.2, y, label, va="center", color=color, fontsize=7)
    if name:
        ax.text(0, .72, name, fontweight="bold", fontsize=8)
    ax.set(xlim=(-.05, 11), ylim=(-.5-.36*len(scorer_rows), 1))
    ax.axis("off")


def radar(ax, profiles, metrics):
    angles = np.linspace(np.pi, 0, len(metrics))
    for radius in [.25, .5, .75, 1.]:
        xx = np.linspace(0, np.pi, 150)
        ax.plot(radius*np.cos(xx), radius*np.sin(xx), color="#D9E0E7", lw=.5)
    for theta, metric in zip(angles, metrics):
        ax.plot([0, np.cos(theta)], [0, np.sin(theta)], color="#D9E0E7", lw=.5)
        ax.text(1.15*np.cos(theta), 1.15*np.sin(theta)-(.04 if theta in (0, np.pi) else 0), metric, ha="center", va="center", fontsize=7)
    for i, profile in enumerate(profiles):
        values = np.asarray(profile["normalized_values"], float)
        mean = values.mean(0)
        fixed = profile["role"] == "initial" and np.all(values == values[0])
        sd = np.zeros(len(metrics)) if fixed else values.std(0, ddof=1)
        color = profile.get("color", "#9AA4B2" if profile["role"] == "initial" else COLORS[(i-1) % len(COLORS)])
        ax.plot(mean*np.cos(angles), mean*np.sin(angles), color=color, lw=1.25, label=profile.get("label", LABELS.get(profile["role"], profile["role"])))
        for theta, m, s in zip(angles, mean, sd):
            ax.plot([(m-s)*np.cos(theta), (m+s)*np.cos(theta)], [(m-s)*np.sin(theta), (m+s)*np.sin(theta)], color=color, lw=.7)
    ax.plot([-1, 1], [0, 0], color="#D9E0E7", lw=.5)
    ax.set(xlim=(-1.5, 1.5), ylim=(-.28, 1.36), aspect="equal")
    ax.axis("off")
    ax.legend(loc="upper center", bbox_to_anchor=(.5, -.02), ncol=2, frameon=False)


def codon_maps(ax, maps):
    for i, item in enumerate(maps):
        pos = item["changed_codon_positions_1_based"]
        role = item["role"]
        color = item.get("color", "#9AA4B2" if role == "initial" else COLORS[(i-1) % len(COLORS)])
        ax.add_patch(Rectangle((0, i-.21), 117, .42, facecolor="#EDF1F5", lw=0))
        for p in pos:
            ax.add_patch(Rectangle((p-1, i-.21), 1, .42, facecolor=color, lw=0))
        ax.text(-3, i, LABELS.get(role, role), ha="right", va="center", fontsize=6)
        ax.text(120, i, f'{len(pos)}/117  {item["five_prime_utr"].removeprefix("5UTR_")} / {item["three_prime_utr"].removeprefix("3UTR_")}', va="center", fontsize=5.5)
    ax.set(xlim=(-4, 225), ylim=(len(maps)-.4, -.7))
    ax.axis("off")


def fig2(output, fmt):
    payload = read_json(DATA / "SupplementaryFigure3_Mauger_U/plot_data.json")
    values = np.asarray(payload["percentiles"])
    fig = plt.figure(figsize=(9.8, 4.4))
    grid = fig.add_gridspec(1, 4, width_ratios=[.65, 7, 1, .9], left=.26, right=.99, top=.9, bottom=.15, wspace=.1)
    member, ax, rho, pval = [fig.add_subplot(grid[0, i]) for i in range(4)]
    labels = [label.replace("TE–Orthrus (updated)", "Orthrus TE").replace("(human ensemble)", "").replace("(human percentile)", "") for label in payload["labels"]]
    cmap = LinearSegmentedColormap.from_list("rank", ["#EFF4F7", "#AFD2D3", "#5BA3AD", "#176B82", "#183E61"])
    image = ax.imshow(values, aspect="auto", cmap=cmap, vmin=0, vmax=100, interpolation="nearest")
    ax.set(xticks=[], yticks=np.arange(len(labels)), yticklabels=labels)
    ax.tick_params(axis="y", length=0, pad=65)
    ax.set_title("Luciferase variant", fontsize=8)
    ax.set_xticks(np.arange(-.5, 39, 1), minor=True)
    ax.set_yticks(np.arange(-.5, len(labels), 1), minor=True)
    ax.grid(which="minor", color="white", lw=.4)
    ax.tick_params(which="minor", bottom=False, left=False)
    statistics = payload["statistics"]
    keys = ["experiment"] + [r["scorer"] for r in statistics]
    for panel in [member, rho, pval]:
        panel.set(xlim=(-.1, 1.1), ylim=(len(labels)-.5, -.5))
        panel.axis("off")
    member.set_title("AR AT", fontsize=7)
    rho.set_title("Spearman ρ", fontsize=7)
    pval.set_title("p-value", fontsize=7)
    groups = payload["combinations"]
    for i, key in enumerate(keys):
        if i and key not in groups:
            for x, group in [(0.15, "All Relevant"), (.8, "all translation")]:
                member.scatter(x, i, marker="s", s=13, facecolor="#586AA1" if key in groups[group] else "white", edgecolor="#9AA4B2", lw=.4)
        if i:
            stat = statistics[i-1]
            rho.text(.5, i, f'{stat["spearman_rho"]:+.2f}', ha="center", va="center")
            pv = stat["p_value"]
            pval.text(.5, i, f"{pv:.3g}" + ("*" if pv < .05 else ""), ha="center", va="center", color="#9A461E" if pv < .05 else INK)
    cax=fig.add_axes([.40,.055,.36,.025])
    fig.colorbar(image, cax=cax, orientation="horizontal", label="Rank percentile", ticks=[0, 25, 50, 75, 100])
    save(fig, output, NAMES[1], fmt)
    return {"heatmap": values.tolist(), "statistics": statistics, "column_ids": payload["ids"]}


def fig3(output, fmt):
    payload = read_json(DATA / "figure_display_coordinates_0105.json")["figure3"]
    fig = plt.figure(figsize=(9.3, 10))
    grid = fig.add_gridspec(4, 2, height_ratios=[.8, 1.7, 1.15, 1.1], hspace=.75, wspace=.5)
    ax = fig.add_subplot(grid[0, :])
    scope(ax, [0, 2, 4], [("MFE", 0, 5, COLORS[0]), ("Orthrus TE", 0, 5, COLORS[1]), ("CSI", 1, 3, COLORS[2]), ("Ribosome load", 1, 3, COLORS[3]), ("Sequence feasibility", 0, 5, "#777777")])
    title(ax, "a", "Design space and scorer inputs")
    ax = fig.add_subplot(grid[1, 0]); radar(ax, payload["radar_profiles"], ["−MFE", "Orthrus TE", "CSI", "Ribosome load"]);title(ax, "b", "Selected candidates")
    ax = fig.add_subplot(grid[1, 1]);codon_maps(ax, payload["maps"]);title(ax, "", "Sequence changes · seed 211")
    series_c = rows(DATA / "Figure3_logged_scores_0100/figures/Figure3_modular_design.c_series.csv")
    series_d = rows(DATA / "Figure3_logged_scores_0100/figures/Figure3_modular_design.d_series.csv")
    generation = [int(row["generation"]) for row in series_c]
    ax = fig.add_subplot(grid[2, 0]); title(ax, "c", "Optimization · seed 211")
    for k, color in zip(["mfe_per_nt", "translation_efficiency_score", "csi_percentile", "ribodecode_translation_log"], COLORS):
        ax.plot(generation, [float(r[k+"_relative_mean"]) for r in series_c], color=color, lw=1, label=LABELS[k])
    ax.set(xlabel="Generation", ylabel="Relative score");ax.legend(frameon=False, ncol=2)
    ax = fig.add_subplot(grid[2, 1]);title(ax, "", "Feasibility and archive HV")
    for key, label, color in [("feasible_fraction", "Passing fraction", "#748A6B"), ("relative_archive_hv", "Archive HV", INK)]:
        ax.plot(generation, [float(r[key]) for r in series_c], label=label, color=color)
    ax.set(xlabel="Generation", ylim=(-.03, 1.04));ax.legend(frameon=False)
    ax = fig.add_subplot(grid[3, :]);title(ax, "d", "Additional recorded scores")
    metrics = [("cai", "CAI"), ("orthrus_te_common_holdout", "Orthrus TE"), ("ribonn_mean_te", "RiboNN"), ("half_life_score", "Orthrus half-life"), ("cds_mfe", "CDS MFE"), ("irna_max_profile_risk", "scanMiR"), ("rinalmo_window_logit_p95", "RiNALMo")]
    for (key, label), color in zip(metrics, COLORS):
        ax.plot(generation, [float(r[key+"_relative_mean"]) for r in series_d], color=color, label=label, lw=1)
    ax.set(xlabel="Generation", ylabel="Relative score", ylim=(-.03, 1.04));ax.legend(frameon=False, ncol=4, loc="upper center", bbox_to_anchor=(.5, -.3))
    save(fig, output, NAMES[2], fmt)
    return {"radar": payload["radar_profiles"], "maps": payload["maps"], "c": series_c, "d": series_d}


def fig4(output, fmt):
    payload = read_json(DATA / "figure_display_coordinates_0109.json")
    fig = plt.figure(figsize=(10, 12))
    grid = fig.add_gridspec(6, 2, height_ratios=[.85, 1.25, 1.1, 1.3, .9, 1.1], hspace=.75, wspace=.4)
    panels = [("4a", "VaxPress and MORSE", [2]), ("4b", "LinearDesign refinement", [0, 2, 4]), ("4c", "RiboDecode refinement", [2])]
    for i, (panel, label, mutable) in enumerate(panels):
        comparison=grid[2*i,:].subgridspec(1,2,wspace=.3)
        ax=fig.add_subplot(comparison[0,0])
        upstream=[("MFE",0,5,COLORS[0]),("CSI",1,3,COLORS[2])] if panel=="4a" else [("CDS MFE",1,3,COLORS[0]),("CAI",1,3,COLORS[5])] if panel=="4b" else [("Ribosome load",1,3,COLORS[3])]
        scope(ax,[2],upstream,name="VaxPress" if panel=="4a" else "LinearDesign" if panel=="4b" else "RiboDecode");title(ax,chr(97+i),label)
        ax=fig.add_subplot(comparison[0,1])
        downstream=upstream if panel=="4a" else [("MFE",0,5,COLORS[0]),("CAI",1,3,COLORS[5]),("Orthrus TE",0,5,COLORS[1]),("Feasibility",0,5,"#777777")] if panel=="4b" else [("Ribosome load",1,3,COLORS[3]),("CSI",1,3,COLORS[2]),("Feasibility",0,5,"#777777")]
        scope(ax,mutable,downstream,name="MORSE")
        ax = fig.add_subplot(grid[2*i+1, 0]);radar(ax, [p for p in payload["radar"] if p["panel"] == panel], ["CAI", "CSI", "−MFE", "Orthrus TE", "RiboDecode", "Half-life"])
        ax = fig.add_subplot(grid[2*i+1, 1])
        if panel == "4a":
            scatter = read_json(DATA / "Figure4a_scatter.json")
            for run in scatter["runs"]:
                ax.scatter([-p["mfe_per_nt"] for p in run["points"]], [p["csi_percentile"] for p in run["points"]], s=4, alpha=.3, color=INK if run["condition"]=="morse" else "#B993C8")
                selected=run["balanced"]
                ax.scatter(-selected["mfe_per_nt"],selected["csi_percentile"],s=28,marker="*" if run["condition"]=="morse" else "D",color=INK if run["condition"]=="morse" else "#B993C8",lw=.3)
            ax.scatter(-scatter["initial"]["mfe_per_nt"],scatter["initial"]["csi_percentile"],marker="x",s=22,color="#9AA4B2")
            ax.set(xlabel="−MFE · kcal/(mol·nt)", ylabel="CSI percentile")
        elif panel == "4b":
            codon_maps(ax, payload["maps"])
        else:
            feasibility=payload["feasibility"]
            for j,row in enumerate(feasibility):
                ax.text(-.3,j,str(row["seed"]),ha="right",va="center",fontsize=7)
                for x,passed in [(0,row["initial_pass"]),(1,row["best_hv_pass"])]:
                    ax.scatter(x,j,marker="o",facecolors="white",edgecolors="#509893" if passed else "#9A461E",s=85,lw=.8)
                    ax.text(x,j,"✓" if passed else "×",ha="center",va="center",color="#509893" if passed else "#9A461E",fontsize=8)
                ax.annotate("",xy=(.8,j),xytext=(.2,j),arrowprops={"arrowstyle":"->","color":"#9AA4B2","lw":.6})
            ax.text(0,-.9,f'{sum(r["initial_pass"] for r in feasibility)}/5 pass',ha="center",fontsize=7)
            ax.text(1,-.9,f'{sum(r["best_hv_pass"] for r in feasibility)}/5 pass',ha="center",fontsize=7)
            ax.set(xlim=(-.6,1.4),ylim=(-1.2,4.7),xticks=[0,1],xticklabels=["Native","Best HV"],yticks=[])
            for spine in ax.spines.values():spine.set_visible(False)
    save(fig, output, NAMES[3], fmt)
    return {"radar": payload["radar"], "maps": payload["maps"], "scatter": read_json(DATA / "Figure4a_scatter.json")}


def fig5(output, fmt):
    payload = read_json(DATA / "figure_display_coordinates_0105.json")["figure5"]
    fig = plt.figure(figsize=(9.3, 8.3)); grid = fig.add_gridspec(3, 2, height_ratios=[.6, 1.2, 1.2], hspace=.7, wspace=.5)
    ax = fig.add_subplot(grid[0, :]);scope(ax, [0, 2], [("RiboNN", 0, 5, COLORS[4]), ("Orthrus TE", 0, 5, COLORS[1]), ("Ribosome load", 1, 3, COLORS[3])]);title(ax, "a", "Design space and predictor inputs")
    ax = fig.add_subplot(grid[1, 0]);radar(ax, payload["radar_profiles"], ["RiboNN", "Orthrus TE", "RiboDecode"]);title(ax, "b", "Specialists and compromise")
    ax = fig.add_subplot(grid[1, 1]);bars = payload["bars"];metrics = list(dict.fromkeys(b["metric"] for b in bars));roles=list(dict.fromkeys(b["role"] for b in bars))
    for j, metric in enumerate(metrics):
        bb=[b for b in bars if b["metric"]==metric];x=np.arange(len(bb))+(j-1)*.24
        ax.bar(x,[b["mean"] for b in bb],width=.23,yerr=[b["sample_sd"] or 0 for b in bb],color=[COLORS[4],COLORS[1],COLORS[3]][j],label=LABELS[metric],error_kw={"lw":.6,"capsize":1.5})
    ax.set(xticks=np.arange(len(roles)),xticklabels=[LABELS[r] for r in roles],ylabel="Empirical percentile",ylim=(0,105));ax.tick_params(axis="x",labelrotation=35);ax.legend(frameon=False)
    ax1=fig.add_subplot(grid[2,0]);ax2=fig.add_subplot(grid[2,1]);title(ax1,"c","Rolling predictor correlations");title(ax2,"","Optimization dynamics")
    for curve in payload["curves"]:
        ax=ax1 if curve["key"].startswith("rho_") else ax2
        mean=np.asarray([np.nan if v is None else v for v in curve["mean"]]);lower=np.asarray([np.nan if v is None else v for v in curve["lower"]]);upper=np.asarray([np.nan if v is None else v for v in curve["upper"]]);x=curve["generations"]
        label=curve["label"].replace(" (held-out)","");ax.plot(x,mean,color=curve["color"],label=label,lw=1);ax.fill_between(x,lower,upper,color=curve["color"],alpha=.1,lw=0)
    for ax in [ax1,ax2]:ax.set(xlabel="Generation");ax.legend(frameon=False,fontsize=5.5)
    ax1.set(ylabel="Spearman ρ",ylim=(-1.02,1.02));ax2.set(ylabel="Relative score / archive HV")
    save(fig,output,NAMES[4],fmt)
    return {k:payload[k] for k in ["radar_profiles","bars","curves","trajectory_normalization"]}


def ed1(output,fmt):
    rr=[r for r in runs() if r["figure"]=="F2"]
    policy=read_json(DATA/"Sequence_feasibility_manufacturing_strategy_v1.json")["parameters"]
    metrics=[("mean_feasibility_pass","Population feasibility",(0,1),None),("mean_feasibility_failed_check_count","Failed check types",None,None),
             ("mean_global_gc_fraction","Global GC",None,(policy["global_gc_min"],policy["global_gc_max"])),("mean_local_gc_min","Minimum local GC",None,(policy["local_gc_min"],)),("mean_local_gc_max","Maximum local GC",None,(policy["local_gc_max"],)),
             ("mean_max_homopolymer_a","Longest A run",None,(policy["homopolymer_max_a"],)),("mean_max_homopolymer_c","Longest C run",None,(policy["homopolymer_max_c"],)),("mean_max_homopolymer_g","Longest G run",None,(policy["homopolymer_max_g"],)),("mean_max_homopolymer_u","Longest U run",None,(policy["homopolymer_max_u"],)),
             ("fraction_failed_direct_repeats","Direct-repeat violations",(0,1),None),("fraction_failed_tandem_repeats","Tandem-repeat violations",(0,1),None),
             ("fraction_failed_restriction_sites","Restriction-site violations",(0,1),None),("fraction_failed_capping_start","Transcript-start violations",(0,1),None),("fraction_failed_inverted_repeats","Inverted-repeat violations",(0,1),None)]
    fig,axes=plt.subplots(4,4,figsize=(11,9));series=[]
    for i,(metric,label,limits,thresholds) in enumerate(metrics):
        ax=axes.flat[i];title(ax,chr(97+i),label)
        for condition,color in [("filtered","#509893"),("logging_only","#B993C8")]:
            subset=[r for r in rr if r["condition"]==condition];support=min(r["generations"] for r in rr);x=np.arange(1,support+1)
            values=np.asarray([[c[metric] for c in r["curves"][:support]] for r in subset]);mean=values.mean(0);sd=values.std(0,ddof=1)
            ax.plot(x,mean,color=color,lw=1,label=condition.replace("_"," "));ax.fill_between(x,mean-sd,mean+sd,color=color,alpha=.12,lw=0)
            for r in subset:ax.plot([c["generation"] for c in r["curves"]],[c[metric] for c in r["curves"]],color=color,alpha=.18,lw=.4)
            series.append({"panel":chr(97+i),"metric":metric,"condition":condition,"common_generations":x.tolist(),"mean":mean.tolist(),"sample_sd":sd.tolist()})
        if limits:ax.set_ylim(*limits)
        for threshold in thresholds or []:ax.axhline(threshold,color="#9AA4B2",lw=.5,ls="--")
        ax.set_xlabel("Generation")
    for ax in list(axes.flat)[14:]:ax.axis("off")
    axes.flat[0].legend(frameon=False,fontsize=5.5);fig.tight_layout(h_pad=2,w_pad=1)
    save(fig,output,NAMES[5],fmt);return series


def ed2(output,fmt):
    rr=rows(DATA/"SupplementaryFigure2_LinearDesign_paired_changes.csv");metrics=list(dict.fromkeys(r["metric"] for r in rr))
    fig,axes=plt.subplots(2,4,figsize=(10,4.8))
    for j,metric in enumerate(metrics):
        subset=[r for r in rr if r["metric"]==metric];a=np.asarray([float(r["cds_only"]) for r in subset]);b=np.asarray([float(r["cds_utrs"]) for r in subset])
        ax=axes[0,j];title(ax,chr(97+j),LABELS.get(metric,"CDS MFE"));ax.bar([0,1],[a.mean(),b.mean()],yerr=[a.std(ddof=1),b.std(ddof=1)],color=["#9AA4B2",COLORS[j]],width=.6,error_kw={"lw":.7,"capsize":2})
        for x,y in zip(a,b):ax.plot([0,1],[x,y],color=INK,alpha=.3,lw=.5,marker=".",ms=3)
        ax.set(xticks=[0,1],xticklabels=["Fixed UTRs","Editable UTRs"]);ax.tick_params(axis="x",labelrotation=25)
        ax.set_ylabel({"mfe_per_nt":"kcal/(mol·nt)","cds_mfe":"kcal/mol","cai_log":"log(CAI)","translation_efficiency_score":"Predicted score"}[metric])
        ax=axes[1,j];title(ax,chr(101+j),"Paired change");d=b-a;ax.axhline(0,color="#B8C3CD",lw=.6);ax.scatter(np.arange(5),d,color=COLORS[j],s=18);ax.set(xticks=np.arange(5),xticklabels=[r["seed"] for r in subset],xlabel="Seed")
    fig.tight_layout(w_pad=2,h_pad=2);save(fig,output,NAMES[6],fmt);return rr


def ed3(output,fmt):
    rr=rows(DATA/"SupplementaryFigure1_RiboDecode_paired_changes.csv");fig,axes=plt.subplots(1,2,figsize=(7,3.2))
    for j,key in enumerate(["csi_percentile","ribodecode_translation_log"]):
        ax=axes[j];title(ax,chr(97+j),"Codon stability" if j==0 else "Predicted ribosome load")
        for i,row in enumerate(r for r in rr if r["metric"]==key):ax.plot([0,1],[float(row["initial"]),float(row["selected"])],marker="o",ms=3,lw=.8,color=COLORS[i],label=f'Seed {row["seed"]}')
        ax.set(xticks=[0,1],xticklabels=["Native start","MORSE continuation"],ylabel="CSI percentile" if j==0 else "Native log score z")
    axes[1].legend(frameon=False,bbox_to_anchor=(1,1));fig.tight_layout(w_pad=3);save(fig,output,NAMES[7],fmt);return rr


def build(output,fmt="both"):
    from .verification import verify
    from .analysis import tables
    verify();style();output=Path(output).resolve();output.mkdir(parents=True,exist_ok=True)
    # The concept-only schematic is itself the explicitly versioned vector source.
    shutil.copyfile(HERE/"reference/Figure1_design_workflow.pdf",output/(NAMES[0]+".pdf"))
    if fmt in {"png","both"}:
        # The concept-only vector has a versioned 144-dpi companion. Reusing it
        # avoids loading a native PDF rasterizer into the numerical plotting
        # process; no numerical observations are rasterized or substituted.
        shutil.copyfile(HERE/"reference/diagrams/Figure1_design_workflow.png",output/(NAMES[0]+".png"))
    payload={"Figure2":fig2(output,fmt),"Figure3":fig3(output,fmt),"Figure4":fig4(output,fmt),"Figure5":fig5(output,fmt),
             "ExtendedDataFigure1":ed1(output,fmt),"ExtendedDataFigure2":ed2(output,fmt),"ExtendedDataFigure3":ed3(output,fmt)}
    required=[output/(name+"."+extension) for name in NAMES for extension in (["pdf","png"] if fmt=="both" else [fmt])]
    missing=[path.name for path in required if not path.is_file() or path.stat().st_size==0]
    if missing:raise RuntimeError(f"Figure generation is incomplete: {missing}")
    if fmt in {"png","both"} and sha256(output/(NAMES[0]+".png"))!=sha256(HERE/"reference/diagrams/Figure1_design_workflow.png"):
        raise RuntimeError("Static diagram preview differs from its reviewed source")
    write_json(output/"plotted_observations.json",payload)
    table_files=tables(output/"tables")
    result={"status":"pass","figures":NAMES,"layout":"Reconstructed scientific layouts; exact published PDFs are archived separately in paper/reference.",
            "numerical_source":"Frozen individual observations, current panel-specific normalizations and selected IDs; no old numerical PDF is used as an input.",
            "static_vector_source":"Figure1_design_workflow.pdf is a concept-only diagram.","tables":table_files,
            "outputs":{p.name:sha256(p) for p in sorted(output.glob("*.pdf"))},"plotted_observations_sha256":sha256(output/"plotted_observations.json"),
            "new_model_inference":False,"optimization_rerun":False,"pixel_identical_claimed":False}
    write_json(output/"figure_receipt.json",result);return result
