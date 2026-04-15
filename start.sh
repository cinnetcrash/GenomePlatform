#!/usr/bin/env bash
# GenomePlatform başlatıcı
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BACKEND="$SCRIPT_DIR/backend"

# API anahtarını kontrol et
if [[ -z "${ANTHROPIC_API_KEY:-}" ]]; then
  echo "⚠️  UYARI: ANTHROPIC_API_KEY ayarlanmamış."
  echo "   AI yorumu devre dışı kalacak."
  echo "   Ayarlamak için: export ANTHROPIC_API_KEY='sk-ant-...'"
  echo ""
fi

# Veri dizinlerini oluştur
mkdir -p "$SCRIPT_DIR/data/"{uploads,results,logs}

# Bağımlılıkları kur (yoksa)
pip install -q -r "$SCRIPT_DIR/requirements.txt"

echo "🚀 GenomePlatform başlatılıyor..."
echo "   Adres: http://localhost:8000"
echo "   Durdurmak için: Ctrl+C"
echo ""

cd "$BACKEND"
python -m uvicorn main:app \
  --host 0.0.0.0 \
  --port 8000 \
  --log-level info \
  --reload
