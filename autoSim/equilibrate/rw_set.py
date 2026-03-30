import sys

def rewrite_settings(filename,):
    with open(filename,'r') as f:
        lines = [line.split() for line in f]

    for line in lines:
        try: 
            if line[0] == '#': continue
        except IndexError: 
            continue
        for i, character in enumerate(line[1:]):
            try:
                float(character)
            except ValueError:
                del line[i+1]
                break

    with open(f'system.in.settings','w') as f:
        for line in lines:
            for _, l in enumerate(line):
                if _ == 0: continue
                try: 
                    float(line[_])
                except ValueError:
                    line[_] = f'#{line[_]}'
                    break
                except IndexError:
                    pass

            f.write(" ".join(map(str, line)) + '\n')
#####

if __name__ == '__main__':
    rewrite_settings('system.in.settings')

