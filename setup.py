from setuptools import setup, find_packages

setup(
    name='QuantumMST',
    version='0.1.0',
    packages=find_packages(),
    install_requires=[
        'numpy',
        'scipy',
        'matplotlib',
        'pandas',
        'qiskit',
        'qiskit-aer',
        'qiskit-algorithms',
        'qiskit-ibm-runtime',
    ],
)