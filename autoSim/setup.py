import setuptools

setuptools.setup(
        name='autoSim',
        version='0.1.0',
        description='Python package to automate MDworkflow for small molecule simulations',
        url='https://github.com/Corcellilab/Phosphates/tree/main/autoSim',
        author='Noah Vasconez',
        author_email='nvascone@nd.edu',
        license='MIT',
        install_requires=[
            'MDAnalysis>=2.10.0',
            'moltemplate>=2.22.3',
            'networkx>=3.4.2',
            'pandas>=2.3.0',
            'rdkit>=2025.9.3',
        ],
        classifiers=[
            'Development Status :: 1 - Planning',
            'Intedend Audience :: Science/Research',
            'Programming Language :: Python :: 3.12.12',
        ],
        packages=setuptools.find_packages(),
        package_dir={'autoSim': 'autoSim'},
        package_data={'autoSim': [
                'parameterize/master_gcrt',
                'parameterize/mol_templates/*.lt',
                'equilibrate/lmp_ins/in.*',
                'prod/lmp_ins/in.*'
            ]
        },
    )

                        

