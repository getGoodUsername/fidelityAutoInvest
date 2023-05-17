import math
import subprocess

# when buying assets make limit order that will always buy at 0.3% lower than market ask price
args = ["./a.out", "345000",  "0.65, 0.06, 0, 0, 0, 0.205, 0.085",  "23672.72, 15918.04, 2222, 12222, 24222, 144047.19, 11317.59"]
result = subprocess.run(args, capture_output=True, text=True);
print(result.stdout, end = "")
print(result.stderr, end = "")

# print(math.inf == float('inf'));
# ' '.join(map(str, [1, 2, 55, 321]))
# list(map(float, "1, 2, 3, 4".split(",")))