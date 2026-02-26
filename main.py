import subprocess
import sys
import os

def run(workflow_name: str):
    path = os.path.join("workflows", workflow_name)

    env = os.environ.copy()
    env["PYTHONPATH"] = os.getcwd()

    print(f"inicializando {workflow_name}")

    subprocess.run(
        [sys.executable, path],
        check=True,
        env=env,
        capture_output=False
    )


if __name__ == '__main__':

    run("workflow_optimizer.py")
