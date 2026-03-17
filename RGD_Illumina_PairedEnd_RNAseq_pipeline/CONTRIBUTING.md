# Contributing

Thank you for your interest in contributing to this pipeline!

## Reporting Issues

Please open a GitHub Issue and include:
- The script name and version (check the script header comments)
- Your HPC environment (OS, SLURM version, module versions)
- The relevant section of the log file
- The full command or sbatch invocation used

## Suggesting Changes

Open a GitHub Issue to discuss before submitting a pull request for significant changes.
For small bug fixes, a pull request with a clear description is welcome directly.

## Pull Request Guidelines

1. Fork the repository and create a branch from `main`
2. Test your changes on at least one real BioProject before submitting
3. Update `CHANGELOG.md` under `[Unreleased]` describing your change
4. Update `README.md` if you change script behavior, inputs, or outputs
5. Keep SLURM directives (`#SBATCH`) generic — users will need to update
   account names, partitions, and email addresses for their own cluster

## Code Style

- Use `set -euo pipefail` at the top of all bash scripts
- Validate required variables and files early, with clear error messages
- Use the shared `log()` function from `lib_v10.sh` for timestamped output
- Avoid hardcoded absolute paths — use the USER CONFIGURATION block
- Comment non-obvious logic, especially around SLURM dependency chains
