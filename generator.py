import glob, os

pages = [
    {
        'id': 'home',
        'title': 'Roman Gerasimov',
        'menu': 'Home',
    },
    {
        'id': 'globs',
        'title': 'Globular clusters',
        'menu': 'Research/Globular clusters',
    },
    {
        'id': 'popIII',
        'title': 'First stars',
        'menu': 'Research/First stars',
    },
    {
        'id': 'qg',
        'title': 'Physics beyond the Standard Model',
        'menu': 'Research/Quantum gravity',
    },
    {
        'id': 'models',
        'title': 'Stellar models',
        'menu': 'Models',
    },
    {
        'id': 'eheu',
        'title': 'Standard solar abundances',
        'menu': 'Standard chemistry',
    },
    {
        'id': 'cv',
        'title': 'Curriculum vitae',
        'menu': 'Curriculum vitae',
    },
]


# Clear existing pages
existing = glob.glob('./*.html')
for page in existing:
    os.remove(page)

# Load template
f = open('template.tmp', 'r')
template = f.read()
f.close()

# Build menu
included = []
menu = '<ul>'
for page in pages:
    if page['menu'].find('/') == -1:
        title = page['menu']
        url = './{}.html'.format(page['id'])
        if page['id'] == 'home':
            url = './'
        menu += '<li><a href="{}">{}</a></li>'.format(url, title)
    else:
        top = page['menu'].split('/')[0]
        if top in included:
            continue
        included += [top]
        menu += '<li><span class="opener">{}</span><ul>'.format(top)
        for subpage in pages:
            if subpage['menu'].find('/') != -1 and subpage['menu'].split('/')[0] == top:
                menu += '<li><a href="./{}.html">{}</a></li>'.format(subpage['id'], subpage['menu'].split('/')[1])
        menu += '</ul>'
menu += '</ul>'

# Generate new pages
for page in pages:
    content = template + ''
    content = content.replace('[TITLE]', page['title'])
    if page['id'] == 'home':
        filename = 'index.html'
    else:
        filename = '{}.html'.format(page['id'])
    body = './pages/{}.html'.format(page['id'])
    if os.path.isfile(body):
        f = open(body, 'r')
        body = f.read()
        f.close()
    else:
        body = '<section><header class="main"><h1>Coming soon!</h1></header>'
    content = content.replace('[CONTENT]', body)
    content = content.replace('[MENU]', menu)
    f = open(filename, 'w')
    f.write(content)
    f.close()

print('All done!')
    
