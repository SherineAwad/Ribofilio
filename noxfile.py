import nox

@nox.session
def tests(session):
    session.install("-r", "requirements.txt")
    session.run("pytest", "/home/ec2-user/ribofilio/tests/tests.py")

@nox.session
def lint(session):
    session.install("flake8")
    session.run("flake8", "src/ribofilio.py")

