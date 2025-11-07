import json

def write_to_json(data, output_path):
  with open(output_path, 'w') as f:
    json.dump(data, f)