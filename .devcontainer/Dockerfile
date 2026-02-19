FROM mcr.microsoft.com/devcontainers/go:2-1.25-bookworm
RUN apt-get update \
    && apt-get install -y vim \
    && apt-get install -y build-essential procps curl file git \
    && apt-get clean \
    && /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)" \
    && echo >> /home/vscode/.bashrc \
    && echo 'eval "$(/home/linuxbrew/.linuxbrew/bin/brew shellenv)"' >> /home/vscode/.bashrc \
    && eval "$(/home/linuxbrew/.linuxbrew/bin/brew shellenv)" \
    && /home/linuxbrew/.linuxbrew/bin/brew install codex

