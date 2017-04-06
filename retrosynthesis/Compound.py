import numpy as np
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import rdDepictor
from IPython.display import Image

class Compound:
    """Coumpound class for Retrosynthesis"""
    def __init__(self, mol):
        self.mol = mol
        self.comments = []
        self.heads = []
        self.atom_block = []
        self.bond_block = []
        self.tails = []
        self.n_atoms = 0
        self.n_bonds = 0
        
        # 複数の原子（０スタート）を引数とし、結合のID（１スタート）を返す。結合がない場合は０を返す。
        self.matrix = np.zeros([mol.GetNumAtoms(), mol.GetNumAtoms()])
        self.visible_bonds = []
        self.atom_symbols = []
        self.InputMolBlock(Chem.MolToMolBlock(mol))
        self.parents = []
        self.children = []
    
    def InputMolBlock(self, molblock):
        for i, line in enumerate(molblock.split("\n")):
            if i < 4:
                #print("heads", line)
                self.heads.append(line)
                if i == 3:
                    n_atoms, n_bonds = [int(n) for n in line.split()[:2]]
            elif i < (4 + n_atoms):
                #print("atoms", line)
                self.atom_block.append(line)
                a = line.split()
                self.atom_symbols.append(a[3])
            elif i < (4 + n_atoms + n_bonds):
                self.bond_block.append(line)
                self.visible_bonds.append(True)
                n1, n2 = line.split()[:2]
                n1, n2 = int(n1) - 1, int(n2) - 1
                self.matrix[n1][n2] = int(len(self.bond_block))
                self.matrix[n2][n1] = int(len(self.bond_block))
            else:
                #print("tail", line)
                #print(line)
                if 'CHG' in line.split(): continue
                self.tails.append(line)
        #self.visible_atoms = [True for n in range(self.n_atoms)]
        self.n_atoms = n_atoms
        self.n_bonds = n_bonds
        
    def GetMolfile(self):
        #print(self.n_bonds, self.visible_bonds)
        out_str = ''
        #out_str += "\n".join(self.heads) + "\n"
        for i, line in enumerate(self.heads):
            if i == 3:
                out_str += "%3d%3d  0  0  0  0  0  0  0  0999 V2000\n" % (self.n_atoms, self.n_bonds)
            else:
                out_str += line + "\n"
        out_str += "\n".join(self.atom_block) + "\n"
        for line, visible in zip(self.bond_block, self.visible_bonds):
            if visible: out_str += line + "\n"
        out_str += "\n".join(self.tails) + "\n"
        return out_str

    def PrintMatrix(self):
        print(self.matrix)

    def FindSeq(self, length):
        #seqs = []
        seen = []
        for i in range(len(self.atom_block)):
            queue = []
            for j in range(len(self.atom_block)):
                if self.matrix[i][j] == 0: continue
                queue.append([i, j])
            while len(queue) > 0:
                #print(queue)
                seq = queue.pop()
                if len(seq) == length:
                    if seq in seen: 
                        continue
                    else:
                        seen.append(seq)
                        yield seq
                    continue
                k = seq[-1]
                for h in range(len(self.atom_block)):
                    if h in seq: continue
                    if self.matrix[k][h] == 0: continue
                    queue.append(seq + [h])

    def FindUniqSeq(self, length):
        seen = []
        for seq in self.FindSeq(length):
            reversed_seq = seq[::-1]
            #print(seq, reversed_seq)
            if seq in seen: continue
            if reversed_seq in seen: continue
            seen.append(seq)
            seen.append(reversed_seq)
            yield seq
            
    def FindRing(self, length):
        #seqs = []
        seen = []
        for i in range(len(self.atom_block)):
            queue = []
            for j in range(len(self.atom_block)):
                if self.matrix[i][j] == 0: continue
                queue.append([i, j])
            while len(queue) > 0:
                #print(queue)
                seq = queue.pop()
                if len(seq) == length + 1:
                    continue
                k = seq[-1]
                for h in range(len(self.atom_block)):
                    if self.matrix[k][h] == 0: continue
                    if h in seq: 
                        if len(seq) == length and seq[0] == h:
                            #print(seq, h)
                            if seq in seen:
                                continue
                            else:
                                seen.append(seq)
                                yield(seq)
                        continue
                    queue.append(seq + [h])

    def FindUniqRing(self, length):
        seen = []
        for seq in self.FindRing(length):
            reversed_seq = seq[::-1]
            #print(seq, reversed_seq)
            if seq in seen: continue
            if reversed_seq in seen: continue
            seen.append(seq)
            seen.append(reversed_seq)
            yield seq

    def DeleteBond(self, i): # 結合のIDを引数とする
        self.visible_bonds[int(i)] = False
        self.n_bonds = sum([1 for n in self.visible_bonds if n == True])

    def AddBond(self, i, j, k): # 二つの原子の間に結合を生成する、あるいは結合の種類を変える。原子の指定は１スタート
        i = int(i)
        j = int(j)
        if self.matrix[i - 1][j - 1] == 0:
            #print(self.bond_block)
            #print(len(self.bond_block))
            self.bond_block.append("%3d%3d%3d  0" % (i, j, k))
            self.visible_bonds.append(True)
            self.matrix[i - 1][j - 1] = int(len(self.bond_block))
            self.matrix[j - 1][i - 1] = int(len(self.bond_block))
            self.n_bonds = sum([1 for n in self.visible_bonds if n == True])
        else:
            self.bond_block[int(self.matrix[i - 1][j - 1] - 1)] = "%3d%3d%3d  0" % (i, j, k)
    
    def SeqBondOrder(self, seq, lst): # 指定した原子列間の結合次数（単結合、二重結合など）が、指定したものであるなら True
        #print(seq, lst)
        for idx in range(len(seq) - 1):
            a1, a2 = seq[idx:idx+2]
            b1 = self.bond_block[int(self.matrix[a1][a2] - 1)].split()
            if float(b1[2]) != float(lst[idx]):
                return False
        return True
    
    def RingBondOrder(self, seq, lst): # 指定した原子列間の結合次数（単結合、二重結合など）が、指定したものであるなら True
        #print(seq, lst)
        for idx in range(len(seq) - 1):
            a1, a2 = seq[idx:idx+2]
            b1 = self.bond_block[int(self.matrix[a1][a2] - 1)].split()
            if float(b1[2]) != float(lst[idx]):
                return False
        a1, a2 = seq[0], seq[-1]
        b1 = self.bond_block[int(self.matrix[a1][a2] - 1)].split()
        if float(b1[2]) != float(lst[-1]):
            return False
        return True