from setuptools import setup
setup(
    name="demux_read_index",
    version="0.2",
    py_modules=['demux_read_index'],
    install_requires=['HTSeq'],
    entry_points='''
       [console_scripts]
       demux_read_index=demux_read_index:main
    '''
)
