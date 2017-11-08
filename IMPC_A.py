import csv
import pandas
import numpy
import sys
import requests
from bs4 import BeautifulSoup
import os
import subprocess
import time
import datetime
import re
import glob
import lxml
from pprint import pprint

directory = 'C:/Users/Branden/Desktop/IMPC/'
# get list of urls to protocol parameter pages
def get_protocol_list(pipeline = "IMPC"):
	if pipeline == "IMPC":
		url = "https://www.mousephenotype.org/impress/procedures/7"
	if pipeline == "JAX":
		url = "https://www.mousephenotype.org/impress/procedures/12"

	# Scrape the HTML at the url
	r = requests.get(url, verify = False)

	# Turn the HTML into a Beautiful Soup object
	soup = BeautifulSoup(r.text, "lxml")

	f = open(directory + 'protocol_list.txt','w+')

		# create a list of links to the protocol pages
	protocol_list = []
	for a in soup.find_all(id=="Display Protocol"):
		if a.has_attr('href'):
			if "protocol" in a['href'] and a['href'] not in protocol_list and "IMPC" in a.text:
				name = a.text.split(" IMPC")[0]
				code = "IMPC" + a.text.split(" IMPC")[1]
				link = a['href']
				protocol_list.append([name, code, link])
				f.write(name + ' ' + code + '\n')
			if "protocol" in a['href'] and a['href'] not in protocol_list and "JAX" in a.text:
				name = a.text.split(" JAX")[0]
				code = "IMPC" + a.text.split(" JAX")[1]
				link = a['href']
				protocol_list.append([name, code, link])
				f.write(name + ' ' + code + '\n')
	return protocol_list

# get procedure group data
def return_procedure_group(procedure_code, number_of_rows, rewrite = False):
		# project directory name
		project_directory = directory + procedure_code

		# create directory if it doesn't exist
		try:
			os.stat(project_directory)
		except:
			os.mkdir(project_directory)

		# output directory name
		output_directory = directory + procedure_code + '/RawData'

		# create directory if it doesn't exist
		try:
			os.stat(output_directory)
		except:
			os.mkdir(output_directory)
		
		# output file name
		output_file = output_directory + '/' + procedure_code + '_data_2.csv'

		# download raw data
		if rewrite or not os.path.isfile(output_file):
			file = open(output_file,'w+')
			print(procedure_code + ": BEGINNING DOWNLOAD OF RAW DATA")
			q = requests.get('http://wwwdev.ebi.ac.uk/mi/impc/dev/solr/experiment/select?q=metadata%3A*&fq=procedure_group:*_' + procedure_code +'&wt=csv&rows=' + str(number_of_rows), verify = False)
			with open(output_file, 'w') as f:
				f.write(q.text)
			print(procedure_code + ": DOWNLOAD COMPLETED")
		else:
			print(procedure_code + ": DATA ALREADY DOWNLOADED")

def DownloadAll(procedure_code, number_of_rows, rewrite = False):
	get_protocol_list(pipeline = "IMPC")
	return_procedure_group(procedure_code, number_of_rows, rewrite)

if __name__ == '__main__':
	DownloadAll("OFD", 1000000000)
