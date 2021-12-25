from Bio import SeqIO
import numpy as np
import matplotlib.pylab as plt


def print_seq(seq):
    s = seq.seq
    m = int(len(s)/10)
    print(seq.description.split('|')[4].strip() + ':')
    for i in range(m):
        if i == m-1:
            print(s[i*10:(i+1)*10], end='\n')
        else:
            print(s[i*10:(i+1)*10], end=" ")
    return 0


def sim(seq1,seq2):
    c = 0
    for i , j in zip(seq1,seq2):
        if i == j:
            c = c+1
    pr = round(c / max(len(seq1),len(seq2)),2) *100
    print(f'Similarity is :%{pr}')
    return 0


def dot_plot(seq1,seq2,win_size):
    m_dot = [[(seq1[i : i+win_size] != seq2[j : j+win_size]) for j in range(len(seq1) - win_size +1)] for i in range(len(seq2) -win_size +1)]
    m_dot = np.array(m_dot)
    m_dot = np.float32(m_dot)
    plt.imshow(m_dot,cmap='gray')
    plt.xlabel("%s (length %i bp)" % ("First sequence", len(seq1)))
    plt.ylabel("%s (length %i bp)" % ("Second sequence", len(seq2)))
    plt.title("Dot matrix using window size %i" % win_size)
    plt.grid()
    plt.show()
    return 0


Seq_1 = ''
Seq_2 = ''
for rec in SeqIO.parse('Sequences/Human_Insulin.fasta','fasta'):
    Seq_1 = rec
for rec in SeqIO.parse('Sequences/Chimpanzee_Insulin.fasta','fasta'):
    Seq_2 = rec
print_seq(Seq_1)
print_seq(Seq_2)
sim(Seq_1.seq,Seq_2.seq)
dot_plot(Seq_1.seq,Seq_2.seq,3)




