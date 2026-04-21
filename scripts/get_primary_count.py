import sys

# Usage: python3 get_primary_count.py spk
# grabs the second line of a flagstat to get primary count reads
with open(sys.argv[1], 'r') as f:
    lines = f.readlines()
print(lines[1].split()[0])
