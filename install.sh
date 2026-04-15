#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
# LycianWay — Interactive Installer
#
# Usage:  bash install.sh
#
# Presents two installation paths:
#   [1] Conda  — installs all conda environments on the host machine
#   [2] Docker — builds and starts the Docker Compose stack
#
# Both paths also set up the .env file if it does not yet exist.
# ─────────────────────────────────────────────────────────────────────────────
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ── Colours ───────────────────────────────────────────────────────────────────
RED='\033[0;31m'; YELLOW='\033[1;33m'; GREEN='\033[0;32m'
CYAN='\033[0;36m'; BOLD='\033[1m'; RESET='\033[0m'

info()    { echo -e "${CYAN}[INFO]${RESET}  $*"; }
ok()      { echo -e "${GREEN}[ OK ]${RESET}  $*"; }
warn()    { echo -e "${YELLOW}[WARN]${RESET}  $*"; }
err()     { echo -e "${RED}[ERR ]${RESET}  $*" >&2; }
header()  { echo -e "\n${BOLD}$*${RESET}"; }


# ── Banner ────────────────────────────────────────────────────────────────────
echo -e "${BOLD}"
echo "  ██╗  ██╗   ██╗ ██████╗██╗ █████╗ ███╗   ██╗██╗    ██╗ █████╗ ██╗   ██╗"
echo "  ██║  ╚██╗ ██╔╝██╔════╝██║██╔══██╗████╗  ██║██║    ██║██╔══██╗╚██╗ ██╔╝"
echo "  ██║   ╚████╔╝ ██║     ██║███████║██╔██╗ ██║██║ █╗ ██║███████║ ╚████╔╝ "
echo "  ██║    ╚██╔╝  ██║     ██║██╔══██║██║╚██╗██║██║███╗██║██╔══██║  ╚██╔╝  "
echo "  ███████╗██║   ╚██████╗██║██║  ██║██║ ╚████║╚███╔███╔╝██║  ██║   ██║   "
echo "  ╚══════╝╚═╝    ╚═════╝╚═╝╚═╝  ╚═╝╚═╝  ╚═══╝ ╚══╝╚══╝ ╚═╝  ╚═╝   ╚═╝  "
echo -e "${RESET}"
echo "  Automated Bacterial Genomic Analysis Platform"
echo "  ─────────────────────────────────────────────"
echo ""


# ── .env setup ────────────────────────────────────────────────────────────────
setup_env() {
    if [[ -f "$SCRIPT_DIR/.env" ]]; then
        ok ".env already exists — skipping."
        return
    fi
    cp "$SCRIPT_DIR/.env.example" "$SCRIPT_DIR/.env"
    warn ".env created from .env.example"
    echo ""
    read -rp "  Enter your Anthropic API key (leave blank to set later): " api_key
    if [[ -n "$api_key" ]]; then
        sed -i "s|sk-ant-YOUR_KEY_HERE|$api_key|" "$SCRIPT_DIR/.env"
        ok "API key saved."
    else
        warn "Remember to set ANTHROPIC_API_KEY in .env before starting."
    fi
}


# ────────────────────────────────────────────────────────────────────────────
#  PATH A — Conda installation
# ────────────────────────────────────────────────────────────────────────────
install_conda() {
    header "=== Conda Installation ==="

    # Detect conda/mamba/micromamba
    MAMBA=""
    for cmd in micromamba mamba conda; do
        if command -v "$cmd" &>/dev/null; then
            MAMBA="$cmd"
            ok "Using: $(command -v $cmd)"
            break
        fi
    done

    if [[ -z "$MAMBA" ]]; then
        err "conda / mamba / micromamba not found."
        echo "  Install Miniforge: https://github.com/conda-forge/miniforge#install"
        exit 1
    fi

    info "System packages (fastp, kraken2) must be installed separately."
    echo "  On Ubuntu/Debian:  sudo apt install fastp kraken2"
    echo "  On RHEL/CentOS:    sudo yum install fastp  (kraken2 via bioconda)"
    echo ""
    read -rp "  Continue with conda environments? [Y/n] " cont
    [[ "${cont:-Y}" =~ ^[Yy]$ ]] || { info "Aborted."; exit 0; }

    # ── Create environments ──────────────────────────────────────────────────
    CHANNELS="-c conda-forge -c bioconda -c defaults"

    create_env() {
        local name="$1"; shift
        if $MAMBA env list 2>/dev/null | grep -qE "^${name}\s"; then
            ok "env '$name' already exists — skipping."
        else
            info "Creating conda env: $name  (packages: $*)"
            $MAMBA create -n "$name" -y $CHANNELS "$@"
            ok "env '$name' created."
        fi
    }

    echo ""
    header "Creating conda environments..."
    echo "  This may take 10–30 minutes depending on your internet connection."
    echo ""

    create_env analiz   "flye>=2.9" "ncbi-amrfinderplus>=3.12" abricate nanoplot bandage
    create_env shovill  shovill
    create_env mlst     mlst
    create_env mobsuite mob_suite
    create_env quast5   "quast>=5.2"
    create_env checkM   checkm2
    create_env checkv   checkv
    create_env bakta    bakta

    # ── Python web-app dependencies ──────────────────────────────────────────
    echo ""
    header "Installing Python web-app dependencies..."
    pip install -q -r "$SCRIPT_DIR/requirements.txt"
    ok "Python packages installed."

    # ── .env ─────────────────────────────────────────────────────────────────
    echo ""
    setup_env

    # ── Done ─────────────────────────────────────────────────────────────────
    echo ""
    ok "Conda installation complete!"
    echo ""
    echo "  To start LycianWay:"
    echo -e "    ${BOLD}bash start.sh${RESET}"
    echo ""
    echo "  Then open: http://localhost:8000"
    echo ""
    warn "Don't forget to set ANTHROPIC_API_KEY in .env for AI interpretation."
}


# ────────────────────────────────────────────────────────────────────────────
#  PATH B — Docker installation
# ────────────────────────────────────────────────────────────────────────────
install_docker() {
    header "=== Docker Installation ==="

    # Check Docker
    if ! command -v docker &>/dev/null; then
        err "Docker not found."
        echo "  Install Docker Engine: https://docs.docker.com/engine/install/"
        exit 1
    fi

    if ! docker compose version &>/dev/null; then
        err "Docker Compose (v2) not found."
        echo "  It usually ships with Docker Desktop."
        echo "  Or install the plugin: https://docs.docker.com/compose/install/"
        exit 1
    fi

    ok "Docker $(docker --version | awk '{print $3}' | tr -d ,)"
    ok "Docker Compose $(docker compose version --short)"

    # ── .env ─────────────────────────────────────────────────────────────────
    echo ""
    setup_env

    # Ask for databases directory
    echo ""
    info "Databases are mounted from your host into the container."
    read -rp "  Enter path to your databases directory [/home/analysis/databases_all]: " db_dir
    db_dir="${db_dir:-/home/analysis/databases_all}"

    if [[ ! -d "$db_dir" ]]; then
        warn "Directory '$db_dir' does not exist — you can create it later."
        warn "Set DATABASES_DIR in .env before running 'docker compose up'."
    else
        # Set in .env
        if grep -q "^DATABASES_DIR=" "$SCRIPT_DIR/.env"; then
            sed -i "s|^DATABASES_DIR=.*|DATABASES_DIR=$db_dir|" "$SCRIPT_DIR/.env"
        else
            echo "DATABASES_DIR=$db_dir" >> "$SCRIPT_DIR/.env"
        fi
        ok "DATABASES_DIR set to: $db_dir"
    fi

    # ── Build ─────────────────────────────────────────────────────────────────
    echo ""
    info "Building Docker image..."
    warn "First build takes 30–60 minutes (conda environments are installed inside)."
    echo ""
    read -rp "  Start the build now? [Y/n] " build_now
    if [[ "${build_now:-Y}" =~ ^[Yy]$ ]]; then
        cd "$SCRIPT_DIR"
        docker compose build
        ok "Docker image built."
    else
        info "Run 'docker compose build' when ready."
        return
    fi

    # ── Start ─────────────────────────────────────────────────────────────────
    echo ""
    read -rp "  Start the platform now? [Y/n] " start_now
    if [[ "${start_now:-Y}" =~ ^[Yy]$ ]]; then
        cd "$SCRIPT_DIR"
        docker compose up -d
        echo ""
        ok "LycianWay is running!"
        echo ""
        echo "  Open:          http://localhost:8000"
        echo "  Stream logs:   docker compose logs -f"
        echo "  Stop:          docker compose down"
    else
        echo ""
        ok "Docker image ready. Start later with:"
        echo -e "    ${BOLD}docker compose up -d${RESET}"
    fi
    echo ""
}


# ── Main menu ─────────────────────────────────────────────────────────────────
echo "  Choose installation method:"
echo ""
echo -e "  ${BOLD}[1] Conda${RESET}  — Install conda environments on this machine"
echo -e "           Best for: servers where conda is already set up"
echo -e "           Databases are accessed from their current locations"
echo ""
echo -e "  ${BOLD}[2] Docker${RESET} — Build a self-contained Docker image"
echo -e "           Best for: clean machines, sharing with others, reproducibility"
echo -e "           Requires: Docker Engine + Docker Compose v2"
echo ""
read -rp "  Your choice [1/2]: " choice

case "$choice" in
    1) install_conda  ;;
    2) install_docker ;;
    *)
        err "Invalid choice: '$choice'"
        echo "  Run the script again and enter 1 or 2."
        exit 1
        ;;
esac
