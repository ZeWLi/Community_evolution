import os,sys


def parse_phobius_out(file):
    with open(file) as f:
        a = f.readlines()[1:]

    out = {}

    for line in a:
        split = line.strip('\n').split('\t')
        gene_cluster = split[0].split(',')[0]

        all_length = []

        try:
            for i in split[1].split(','):
                length = 0
                length += eval(i.split('-')[1]) - eval(i.split('-')[0])
                all_length.append(['intracellular',length])

        except IndexError:
            all_length.append(['intracellular', 0])

        try:
            for i in split[2].split(','):
                length = 0
                length += eval(i.split('-')[1]) - eval(i.split('-')[0])
                all_length.append(['extracellular',length])

        except IndexError:
            all_length.append(['extracellular', 0])

        try:
            for i in split[3].split(','):
                length = 0
                length += eval(i.split('-')[1]) - eval(i.split('-')[0])
                all_length.append(['transmembrane',length])

        except IndexError:
            all_length.append(['transmembrane', 0])

        try:
            for i in split[4].split(','):
                length = 0
                length += eval(i.split('-')[1]) - eval(i.split('-')[0])
                all_length.append(['signal',length])

        except IndexError:
            all_length.append(['signal', 0])

        all_length.sort(key=lambda x:x[1],reverse=True)

        out[gene_cluster] = all_length[0][0]

    return out



def main():
    lis = os.listdir(sys.argv[1])
    pwd = os.getcwd()
    os.chdir(sys.argv[1])
    out = {}

    for i in lis:
        os.chdir(i)
        out[i] = parse_phobius_out('phobius.parse.out')
        os.chdir('..')

    os.chdir(pwd)

    with open(sys.argv[2],'w') as f:
        f.write('cluster\tgene_cluster\tsublocation\n')
        for i in out:
            for j in out[i]:
                f.write(i + '\t' + j + '\t' + out[i][j] + '\n')


if __name__ == '__main__':
    main()
