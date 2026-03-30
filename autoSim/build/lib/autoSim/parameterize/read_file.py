from pathlib import Path

####
def read_file(filename):
    current_dir = Path(__file__).parent

    file_path = current_dir / filename

    with open(file_path, 'r', encoding='utf-8') as f:
        file = [line.split() for line in f]
        return file

