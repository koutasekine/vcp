#!/usr/bin/env bash
set -euo pipefail

need_command() {
    if ! command -v "$1" >/dev/null 2>&1; then
        echo "[ERROR] '$1' is required." >&2
        exit 1
    fi
}

need_command tar
need_command wget

if [ "$#" -ge 1 ]; then
    foldername="$1"
else
    echo "Create a new folder in your home directory."
    printf "Please input a new folder name: "
    read -r foldername
fi

if [ -z "${foldername}" ]; then
    echo "[ERROR] Folder name is empty." >&2
    exit 1
fi

case "${foldername}" in
    */*|*\\*)
        echo "[ERROR] Please input only a folder name, not a path." >&2
        exit 1
        ;;
esac

folderpath="${HOME}/${foldername}"

if [ -e "${folderpath}" ]; then
    echo "[ERROR] '${folderpath}' already exists." >&2
    echo "Please choose a new folder name." >&2
    exit 1
fi

tmpdir="$(mktemp -d)"
cleanup() {
    rm -rf "${tmpdir}"
}
trap cleanup EXIT

echo
echo "######################################################"
echo "Download kv library:"
cd "${tmpdir}"

kv_file="kv-master.tar.gz"
kv_dir="kv-master"
kv_url="https://github.com/mskashi/kv/archive/refs/heads/master.tar.gz"

wget -O "${kv_file}" "${kv_url}"
tar -xzf "${kv_file}"

if [ ! -d "${kv_dir}/kv" ]; then
    echo "[ERROR] kv archive does not contain '${kv_dir}/kv'." >&2
    exit 1
fi

echo "######################################################"
echo "Download VCP library:"

vcp_file="vcp-master.tar.gz"
vcp_dir="vcp-master"
vcp_url="https://github.com/koutasekine/vcp/archive/refs/heads/master.tar.gz"

wget -O "${vcp_file}" "${vcp_url}"
tar -xzf "${vcp_file}"

if [ ! -d "${vcp_dir}/vcp" ]; then
    echo "[ERROR] VCP archive does not contain '${vcp_dir}/vcp'." >&2
    exit 1
fi

mkdir -p "${folderpath}"

cp -R "${kv_dir}/kv" "${folderpath}/"
if [ -d "${kv_dir}/test" ]; then
    cp -R "${kv_dir}/test" "${folderpath}/"
fi
if [ -d "${kv_dir}/example" ]; then
    cp -R "${kv_dir}/example" "${folderpath}/"
fi

cp -R "${vcp_dir}/vcp" "${folderpath}/"
if [ -d "${vcp_dir}/test_matrix" ]; then
    cp -R "${vcp_dir}/test_matrix" "${folderpath}/"
fi
if [ -d "${vcp_dir}/test_PDE" ]; then
    cp -R "${vcp_dir}/test_PDE" "${folderpath}/"
fi
if [ -d "${vcp_dir}/vblas" ]; then
    cp -R "${vcp_dir}/vblas" "${folderpath}/"
fi
if [ -d "${vcp_dir}/docs" ]; then
    cp -R "${vcp_dir}/docs" "${folderpath}/"
fi
if [ -d "${vcp_dir}/tools" ]; then
    cp -R "${vcp_dir}/tools" "${folderpath}/"
fi
if [ -f "${vcp_dir}/README.md" ]; then
    cp "${vcp_dir}/README.md" "${folderpath}/"
fi
if [ -f "${vcp_dir}/LICENSE" ]; then
    cp "${vcp_dir}/LICENSE" "${folderpath}/"
fi

echo
echo "Finish!"
echo "Created: ${folderpath}"
echo
echo "Basic matrix example:"
echo "  cd ${folderpath}/test_matrix"
echo "  g++ -std=c++11 -I.. test_matrix.cpp"
echo
echo "If you use BLAS/LAPACK, check that BLAS respects rounding mode changes."
echo "Example with OpenBLAS:"
echo "  cd ${folderpath}/test_matrix"
echo "  g++ -I.. -std=c++11 -DNDEBUG -DKV_FASTROUND -O3 Check_pdblas_rounding.cpp -llapack -lopenblas -lmpfr -lgmp -fopenmp && ./a.out"
echo "  g++ -I.. -std=c++11 -DNDEBUG -DKV_FASTROUND -O3 Check_pidblas_rounding.cpp -llapack -lopenblas -lmpfr -lgmp -fopenmp && ./a.out"
echo
echo "For Intel MKL or other BLAS/LAPACK settings, see:"
echo "  ${folderpath}/docs/blas_build.md"
