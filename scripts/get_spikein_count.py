import sys

# Usage: python3 get_primary_count.py spk
# grabs the third value from the second line as orthogonal reads
with open(sys.argv[1], 'r') as spike_in_f:
    spike_in_f.readline()  # skip header
    summary_line = spike_in_f.readline()
    summary_fields = summary_line.split('\t')
    scale_factor_raw = int(summary_fields[2])
print(scale_factor_raw)
