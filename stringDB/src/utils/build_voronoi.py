import os
import platform
import subprocess

src = "src/utils/voronoiMain2023.cpp"

if platform.system() == "Windows":
    ext = ".dll"
else:
    ext = ".so"

out = f"src/utils/voronoi_c{ext}"
cxx = os.environ.get("CXX", "g++")

cmd = [cxx, "-O3", "-shared", "-o", out, src]

if platform.system() != "Windows":
    cmd.insert(3, "-fPIC")

print("Running:", " ".join(cmd))
subprocess.check_call(cmd)
print("Built:", out)