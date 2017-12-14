
all = ['like_','em_','gibbs_']

flag = open('./flag')
flag = flag.readline()
print(flag)
if flag == '1\n':

	for each in all:
		f = open('./'+each+'predictedSequence.csv')
		w = open('./'+'final_seq'+each+'.csv','w')

		text = f.readline()
		l_text = text.split(',')
		final = '\n'.join(l_text)

		w.write(final)
		f.close()
		w.close()
else:
	for each in all:
		f = open('./'+each+'predictedSequence.csv')
		w = open('./'+'final_seq'+each+'.csv','w')

		text = f.readline()
		l_text = text.split(',')
		line = []
		count = 0
		for word in l_text:
			if count<24:
				count += 1
				line.append(word)
			else:
				w.write(','.join(line)+'\n')
				line = []
				count = 0

		f.close()
		w.close()
		

		