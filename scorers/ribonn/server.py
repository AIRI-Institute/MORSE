import copy
import json
import sys
import traceback


def handle_request(predictor, request):
    unknown = set(request) - {'request_id', 'sequences', 'sequence_contexts', 'parameters', 'params'}
    if unknown:
        raise ValueError('Unknown request fields: ' + ', '.join(sorted(unknown)))
    sequences, contexts = request.get('sequences'), request.get('sequence_contexts')
    if not isinstance(sequences, list) or not all(isinstance(s, str) for s in sequences):
        raise ValueError('sequences must be a list of strings')
    if not isinstance(contexts, list) or len(sequences) != len(contexts) or not all(isinstance(c, dict) for c in contexts):
        raise ValueError('One explicit sequence context is required per sequence')
    if 'parameters' in request and 'params' in request:
        raise ValueError('Use either parameters or params, not both')
    raw = request.get('parameters', request.get('params', {}))
    if not isinstance(raw, dict):
        raise ValueError('parameters must be a mapping')
    unknown = set(raw) - {'batch_size', 'cpu_threads'}
    if unknown:
        raise ValueError('Unknown RiboNN parameters: ' + ', '.join(sorted(map(str, unknown))))
    parameters = {'batch_size': raw.get('batch_size', predictor.batch_size),
                  'cpu_threads': raw.get('cpu_threads', getattr(predictor, 'cpu_threads', 2))}
    for name, value in parameters.items():
        if type(value) is not int or value < 1:
            raise ValueError(name + ' must be a positive integer')
    configured = copy.copy(predictor)
    configured.batch_size = parameters['batch_size']
    configured.cpu_threads = parameters['cpu_threads']
    return {'request_id': request.get('request_id'), 'parameters': parameters,
            'outputs': configured.predict(sequences, contexts)}


def main():
    from model import Predictor
    predictor = Predictor()
    print(f'RiboNN ready: human, 50 models, 78 targets, {predictor.device}', file=sys.stderr, flush=True)
    for line in sys.stdin:
        if not line.strip():
            continue
        request = {}
        try:
            request = json.loads(line)
            response = handle_request(predictor, request)
        except Exception as exc:
            traceback.print_exc(file=sys.stderr)
            response = {'request_id': request.get('request_id') if isinstance(request, dict) else None,
                        'error': str(exc)}
        print(json.dumps(response, allow_nan=False), flush=True)


if __name__ == '__main__':
    main()

