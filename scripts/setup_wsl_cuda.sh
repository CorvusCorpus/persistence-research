#!/bin/bash
# Set up a dedicated WSL2 distribution for this repository: compiler, GMP/MPFR, Python packages, CUDA compiler.
# Used on the desktop (Windows 11, RTX 3060 Ti) in session 3; see docs/NOTEBOOK.md, task 1 on the desktop.
#
# From Windows (PowerShell), with the Ubuntu 24.04.5 WSL image from releases.ubuntu.com:
#   wsl --install --from-file ubuntu-24.04.5-wsl-amd64.wsl --name persistence --no-launch
#   wsl -d persistence -u root -- bash scripts/setup_wsl_cuda.sh        (from the repository root)
#   wsl --terminate persistence                                          (applies /etc/wsl.conf)
# Afterwards run everything through a login shell, which puts the venv and CUDA on PATH:
#   wsl -d persistence -u root -- bash -lc 'cd /mnt/c/.../persistence-research && make'
#
# The GPU driver is the Windows one (/usr/lib/wsl/lib/libcuda.so). Never install a Linux NVIDIA driver in WSL;
# NVIDIA's wsl-ubuntu repository contains the toolkit only.
set -euo pipefail
CUDA_VER=${CUDA_VER:-13-3}                 # CUDA 13.3 matches the Windows driver 610.88 (CUDA UMD 13.3)
KEYRING_SHA256=eea6cc5f0eaeb99082d054b8c05ae206a378e31e88048df0310d59f651dceed2   # cuda-keyring_1.1-1_all.deb
export DEBIAN_FRONTEND=noninteractive

# Linux tools only on PATH (no Windows PATH entries), systemd as in the stock image
cat > /etc/wsl.conf <<'EOF'
[boot]
systemd=true
[interop]
appendWindowsPath=false
EOF

apt-get update
apt-get install -y build-essential libgmp-dev libmpfr-dev python3-venv python3-dev git xz-utils wget ca-certificates

# Python packages in a venv (Ubuntu 24.04 marks the system Python as externally managed)
python3 -m venv /opt/venv
/opt/venv/bin/pip install --upgrade pip
/opt/venv/bin/pip install python-flint mpmath numpy

# CUDA compiler and profiler from NVIDIA's WSL repository
wget -q -O /tmp/cuda-keyring_1.1-1_all.deb \
    https://developer.download.nvidia.com/compute/cuda/repos/wsl-ubuntu/x86_64/cuda-keyring_1.1-1_all.deb
echo "$KEYRING_SHA256  /tmp/cuda-keyring_1.1-1_all.deb" | sha256sum -c -
dpkg -i /tmp/cuda-keyring_1.1-1_all.deb
apt-get update
apt-get install -y cuda-nvcc-$CUDA_VER cuda-cudart-dev-$CUDA_VER cuda-cuobjdump-$CUDA_VER cuda-nvdisasm-$CUDA_VER \
    cuda-profiler-api-$CUDA_VER cuda-nsight-compute-$CUDA_VER

cat > /etc/profile.d/persistence.sh <<EOF
export PATH=/opt/venv/bin:/usr/local/cuda-${CUDA_VER/-/.}/bin:\$PATH
EOF

# record the environment
. /etc/profile.d/persistence.sh
gcc --version | head -1
python3 -c 'import sys, flint, mpmath, numpy; print("python", sys.version.split()[0], "python-flint", flint.__version__, "mpmath", mpmath.__version__, "numpy", numpy.__version__)'
nvcc --version | tail -2
