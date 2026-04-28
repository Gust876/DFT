import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))


def parse_args():
    parser = argparse.ArgumentParser(
        description="Automação de cálculos DFT — otimização geométrica e frequências vibracionais"
    )
    parser.add_argument(
        "workflow",
        choices=["optimizer", "frequencies"],
        help=(
            "'optimizer' para otimização geométrica, "
            "'frequencies' para frequências vibracionais"
        )
    )
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()

    if args.workflow == "optimizer":
        import workflows.workflow_optimizer
    elif args.workflow == "frequencies":
        import workflows.workflow_frequencies
