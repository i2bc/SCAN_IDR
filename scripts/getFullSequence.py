import os, re
import urllib.parse
import urllib.request
import requests
import sys
import argparse

def main(FLAGS):
    
    url = 'https://www.uniprot.org/uploadlists/'

    #print(FLAGS.input)

    for ID in FLAGS.input.split(','):
        print('Getting',ID)
        params = {
        'from':'ACC+ID',
        'to':'ACC',
        'format':'fasta',
        'query': ID,
        'columns' : 'id,entry name,reviewed,protein names,genes,organism,length,database(RefSeq),go(molecular function)',
                }
        fout = open(os.path.join(FLAGS.output_directory, ID+'.fasta'),"w")
        """
            'columns' : 'id,entry name,reviewed,protein names,genes,organism,length,database(RefSeq),go(molecular function)',
            'Content_Type':'form-data'
        """
        url = f'https://rest.uniprot.org/uniprotkb/{ID}?format=fasta'
        fasta = requests.get(url).text
        fout.write(fasta)
        fout.close()
        """
        data = urllib.parse.urlencode(params)
        data = data.encode('utf-8')
        req = urllib.request.Request(url, data)
        with urllib.request.urlopen(req) as f:
           response = f.read()
        """

        
    
if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument(
    '-i',
    '--input',
    type=str,
    required=True,
    help='Retrieve fasta from comma-separated IDs. Example ESA1_YEAST,EPL1_YEAST'
  )
  parser.add_argument(
    '-d',
    '--output_directory',
    type=str,
    default='.',
    help='Directory where to stode the fasta output file'
  )
  FLAGS = parser.parse_args()
  main(FLAGS)





