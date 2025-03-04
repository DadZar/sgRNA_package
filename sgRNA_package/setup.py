from setuptools import setup, find_packages

setup(
    name='sgRNA_package',
    version='1.0',
    packages=find_packages(),
    install_requires=[
        'biopython==1.83',
        'numpy==1.24.4',
        'pandas==2.0.3',
        'tensorflow==2.10.1',
        'joblib==1.4.2',
        'reportlab==3.6.12',
        'xgboost==2.1.4',
        'scikit-learn==1.3.2'
    ],
    entry_points={
        'console_scripts': [
            'sgRNA-run=sgRNA.main:main',
        ]
    },
    description='Paquete para el diseño y predicción de eficiencia de sgRNAs para CRISPR-Cas9.',
    author='Juan David Aguilar Rico & Mariana Valentina Jaimes',
    author_email='davidr.jda@gmail.com',
    url='https://github.com/DadZar/sgRNA_package',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.10',
)
