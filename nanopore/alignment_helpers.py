def get_max_from_file(file):
    max_values = {}
    with open(file) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 4:
                continue
            gene = parts[0]
            max_value = float(parts[3])
            max_values[gene] = (max_value)
    return max_values
