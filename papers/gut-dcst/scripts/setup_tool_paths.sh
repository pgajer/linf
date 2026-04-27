#!/bin/zsh

# Ensure Homebrew and TeX tools are available even in non-login shell sessions.
for tool_dir in /opt/homebrew/bin /usr/local/bin /Library/TeX/texbin; do
  [[ -d "$tool_dir" ]] || continue
  case ":$PATH:" in
    *":$tool_dir:"*) ;;
    *) PATH="$tool_dir:$PATH" ;;
  esac
done

export PATH
