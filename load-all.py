# increase timeout seconds
def loadurl(url, timeout=5):
    temp_name = tmp_filename() + '.' + os.path.splitext(url)[1][1:]

    from urllib.request import urlopen
    content = urlopen(url, timeout=timeout)
    with open(temp_name, 'wb') as f:
        f.write(content.read())
    sage.repl.load.load(temp_name, globals())
    os.unlink(temp_name)

def load_all(isepg=True, timeout=5, load_func='load', local=False):
    # load_func can be 'loadurl' (default) or 'load'
    # local can be False, True, or 'your_path_to_libraries' (e.g., '~/Download/')
    if load_func == 'load':
        func = load
    if load_func == 'loadurl':
        func = lambda url: loadurl(url, timeout=timeout)

    sage_ver = float(sage.misc.banner.SAGE_VERSION)
    
    if isepg:
        ### temporary setting
        ### hope the version become more consistent in the future
        if sage_ver >= 8.9:
            URL = 'https://raw.githubusercontent.com/hunnellm/isepg/refs/heads/main/'
        else:
            URL = 'https://raw.githubusercontent.com/hunnellm/isepg/refs/heads/main/'

            
        files = ['isepg.py','symplectic_eigenvalues.sage']
        for f in files:
            print("Loading %s..."%f);
            func(URL+f) 
