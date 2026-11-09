import json
import sys
import traceback


def handle_request(predictor, request):
    unknown = set(request) - {'request_id', 'sequences', 'sequence_contexts'}
    if unknown:
        raise ValueError('Unknown request fields: ' + ', '.join(sorted(unknown)))
    sequences, contexts = request.get('sequences'), request.get('sequence_contexts')
    if not isinstance(sequences, list) or not all(isinstance(s, str) for s in sequences):
        raise ValueError('sequences must be a list of strings')
    if not isinstance(contexts, list) or len(sequences) != len(contexts) or not all(isinstance(c, dict) for c in contexts):
        raise ValueError('One explicit sequence context is required per sequence')
    return {'request_id': request.get('request_id'), 'outputs': predictor.predict(sequences, contexts)}


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

