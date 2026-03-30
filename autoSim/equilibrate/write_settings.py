import json

settings = {
        'defrost':
            {
                'run': 400000,
                'temp': 300,
                'press': 1,
            },
        'npt':
            {
                'run': 6000000,
                'temp': 300,
                'press': 1,
            },
        'nvt':
            {
                'run': 5000000,
                'temp': 300,
                'press': 1,
            },
        'nve':
            {
                'run': 5000000,
                'temp': 300,
                'press': 1,
            },
    }

with open('settings.json','w') as j:
    json.dump(settings, j, indent=4)

with open('settings.json', 'r') as j:
    data = json.load(j)

print(data['defrost']['run'])

