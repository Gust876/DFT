FROM python:3.10-slim

# Dependências do sistema
RUN apt-get update && apt-get install -y \
    build-essential \
    libopenblas-dev \
    liblapack-dev \
    && rm -rf /var/lib/apt/lists/*

# Diretório de trabalho
WORKDIR /app

RUN pip install pyscf numpy

# Código
COPY . .

CMD ["python", "main.py"]
