import sys
import copy
import numpy as np
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import rdDepictor
from IPython.display import Image
#import Retrosynthesis.Compound
from retrosynthesis.Compound import Compound
class Retrosynthesis:
    """This class uses the Compound class"""
    def __init__(self):
        self.cpds = []
        self.inchis = []
        
    def InputMol(self, mol):
        c = Compound(mol)
        c.comments.append('__init__')
        self.cpds.append(c)
        self.inchis.append(Chem.MolToInchi(mol))
        return True
            
    def GetMolBlock(self, i):
        return Chem.MolToMolBlock(self.cpds[i].mol)
    
    def GetComments(self, i):
        return self.cpds[i].comments
    
    def DrawCpd(self, i, imagefile):
        rdDepictor.Compute2DCoords(self.cpds[i].mol)
        Draw.MolToFile(self.cpds[i].mol, imagefile)
        return True

    def DehydroBondFormation(self, i): # 脱水素的結合生成反応 A-B <= A + B
        for idx, bond in enumerate(self.cpds[i].bond_block):
            new_cpd = copy.deepcopy(self.cpds[i])
            try:
                new_cpd.DeleteBond(idx)
                new_mol = Chem.MolFromMolBlock(new_cpd.GetMolfile())
                new_inchi = Chem.MolToInchi(new_mol)
                new_cpd = Compound.Compound(new_mol)
                new_cpd.comments.append("%d DehydroBondFormation %d" % (i, idx))
            except:
                print(sys.exc_info())
                continue
            if new_inchi in self.inchis: 
                continue
            self.inchis.append(new_inchi)
            self.cpds.append(new_cpd)
            self.cpds[i].children.append(len(self.cpds) - 1)
            self.cpds[len(self.cpds) - 1].parents.append(i)
        return True
    
    def HydroBondDigestion(self, i): # 水素化的結合切断反応 A + B <= A-B
        for idx1, atom1, in enumerate(self.cpds[i].atom_block):
            for idx2, atom2, in enumerate(self.cpds[i].atom_block):
                if idx1 >= idx2:
                    continue
                if self.cpds[i].matrix[idx1][idx2] != 0:
                    continue
                new_cpd = copy.deepcopy(self.cpds[i])
                try:
                    new_cpd.AddBond(idx1, idx2, 1)
                    new_mol = Chem.MolFromMolBlock(new_cpd.GetMolfile())
                    new_inchi = Chem.MolToInchi(new_mol)
                    new_cpd = Compound.Compound(new_mol)
                    new_cpd.comments.append("%d HydroBondDigestion %d %d" % (i, idx1, idx2))
                except:
                    print(sys.exc_info())
                    continue
                if new_inchi in self.inchis: 
                    continue
                self.inchis.append(new_inchi)
                self.cpds.append(new_cpd)
                self.cpds[i].children.append(len(self.cpds) - 1)
                self.cpds[len(self.cpds) - 1].parents.append(i)
        return True
    
    def UnsaturateTransfer2(self, i): # 不飽和結合転位 A-B=C <= A=B-C
        for seq in self.cpds[i].FindSeq(3):
            if self.cpds[i].SeqBondOrder(seq, [1, 2]):
                pass
            else:
                continue
            new_cpd = copy.deepcopy(self.cpds[i])
            try:
                new_cpd.AddBond(seq[0] + 1, seq[1] + 1, 2)
                new_cpd.AddBond(seq[1] + 1, seq[2] + 1, 1)
                new_mol = Chem.MolFromMolBlock(new_cpd.GetMolfile())
                new_inchi = Chem.MolToInchi(new_mol)
                new_cpd = Compound.Compound(new_mol)
                new_cpd.comments.append("%d UnsaturateTransfer2 %d %d %d" % (i, seq[0], seq[1], seq[2]))
            except:
                print(sys.exc_info())
                continue
            if new_inchi in self.inchis: 
                continue
            self.inchis.append(new_inchi)
            self.cpds.append(new_cpd)
            self.cpds[i].children.append(len(self.cpds) - 1)
            self.cpds[len(self.cpds) - 1].parents.append(i)
        return True

    def UnsaturateTransfer4(self, i): # 不飽和結合転位 A-B=C-D=E <= A=B-C=D-E
        for seq in self.cpds[i].FindSeq(5):
            if self.cpds[i].SeqBondOrder(seq, [1, 2, 1, 2]):
                pass
            else:
                print(sys.exc_info())
                continue
            new_cpd = copy.deepcopy(self.cpds[i])
            try:
                new_cpd.AddBond(seq[0] + 1, seq[1] + 1, 2)
                new_cpd.AddBond(seq[1] + 1, seq[2] + 1, 1)
                new_cpd.AddBond(seq[2] + 1, seq[3] + 1, 2)
                new_cpd.AddBond(seq[3] + 1, seq[4] + 1, 1)
                new_mol = Chem.MolFromMolBlock(new_cpd.GetMolfile())
                new_inchi = Chem.MolToInchi(new_mol)
                new_cpd = Compound.Compound(new_mol)
                new_cpd.comments.append(
                    "%d UnsaturateTransfer4 %d %d %d %d %d" % (i, seq[0], seq[1], seq[2], seq[3], seq[4]))
            except:
                print(sys.exc_info())
                continue
            if new_inchi in self.inchis: 
                continue
            self.inchis.append(new_inchi)
            self.cpds.append(new_cpd)
            self.cpds[i].children.append(len(self.cpds) - 1)
            self.cpds[len(self.cpds) - 1].parents.append(i)
        return True

    def UnsaturateTransfer6(self, i): # 不飽和結合転位 A-B=C-D=E-F=G <= A=B-C=D-E=F-G
        for seq in self.cpds[i].FindSeq(7):
            if self.cpds[i].SeqBondOrder(seq, [1, 2, 1, 2, 1, 2]):
                pass
            else:
                continue
            new_cpd = copy.deepcopy(self.cpds[i])
            try:
                new_cpd.AddBond(seq[0] + 1, seq[1] + 1, 2)
                new_cpd.AddBond(seq[1] + 1, seq[2] + 1, 1)
                new_cpd.AddBond(seq[2] + 1, seq[3] + 1, 2)
                new_cpd.AddBond(seq[3] + 1, seq[4] + 1, 1)
                new_cpd.AddBond(seq[4] + 1, seq[5] + 1, 2)
                new_cpd.AddBond(seq[5] + 1, seq[6] + 1, 1)
                new_mol = Chem.MolFromMolBlock(new_cpd.GetMolfile())
                new_inchi = Chem.MolToInchi(new_mol)
                new_cpd = Compound.Compound(new_mol)
                new_cpd.comments.append(
                    "%d UnsaturateTransfer6 %d %d %d %d %d %d %d" % (i, seq[0], seq[1], seq[2], seq[3], seq[4],
                                                                    seq[5], seq[6]))
            except:
                print(sys.exc_info())
                continue
            if new_inchi in self.inchis: 
                continue
            self.inchis.append(new_inchi)
            self.cpds.append(new_cpd)
            self.cpds[i].children.append(len(self.cpds) - 1)
            self.cpds[len(self.cpds) - 1].parents.append(i)
        return True

    def AromaticRingRotate6(self, i): # ベンゼン環中の不飽和結合表記の移動（反応ではない）
        # これを単独で動かしても表記は変わらないので、他の不飽和結合転位と組み合わせる必要性
        for seq in self.cpds[i].FindUniqRing(6):
            if self.cpds[i].RingBondOrder(seq, [1, 2, 1, 2, 1, 2]):
                new_cpd = copy.deepcopy(self.cpds[i])
                #print("self.cpds[i].SeqBondOrder(seq, [1, 2, 1, 2, 1, 2])")
                try:
                    new_cpd.AddBond(seq[0] + 1, seq[1] + 1, 2)
                    new_cpd.AddBond(seq[1] + 1, seq[2] + 1, 1)
                    new_cpd.AddBond(seq[2] + 1, seq[3] + 1, 2)
                    new_cpd.AddBond(seq[3] + 1, seq[4] + 1, 1)
                    new_cpd.AddBond(seq[4] + 1, seq[5] + 1, 2)
                    new_cpd.AddBond(seq[5] + 1, seq[0] + 1, 1)
                    new_mol = Chem.MolFromMolBlock(new_cpd.GetMolfile())
                    new_inchi = Chem.MolToInchi(new_mol)
                    new_cpd = Compound.Compound(new_mol)
                    new_cpd.comments.append(
                        "%d AromaticRingRotate6 %d %d %d %d %d %d" % (i, seq[0], seq[1], seq[2], seq[3], seq[4],
                                                                        seq[5]))
                except:
                    print(sys.exc_info())
                    continue
                if new_inchi in self.inchis: 
                    print("new_inchi in self.inchis")
                    #continue
                else:
                    print("new_inchi NOT in self.inchis")
                    print(len(self.cpds), len(self.inchis), 9999999)
                self.inchis.append(new_inchi)
                self.cpds.append(new_cpd)
                self.cpds[i].children.append(len(self.cpds) - 1)
                self.cpds[len(self.cpds) - 1].parents.append(i)
            elif self.cpds[i].RingBondOrder(seq, [2, 1, 2, 1, 2, 1]):
                print("hoge")
                new_cpd = copy.deepcopy(self.cpds[i])
                #print("self.cpds[i].SeqBondOrder(seq, [2, 1, 2, 1, 2, 1]):")
                try:
                    new_cpd.AddBond(seq[0] + 1, seq[1] + 1, 1)
                    new_cpd.AddBond(seq[1] + 1, seq[2] + 1, 2)
                    new_cpd.AddBond(seq[2] + 1, seq[3] + 1, 1)
                    new_cpd.AddBond(seq[3] + 1, seq[4] + 1, 2)
                    new_cpd.AddBond(seq[4] + 1, seq[5] + 1, 1)
                    new_cpd.AddBond(seq[5] + 1, seq[0] + 1, 2)
                    #print(new_cpd.GetMolFile())
                    new_mol = Chem.MolFromMolBlock(new_cpd.GetMolfile())
                    new_inchi = Chem.MolToInchi(new_mol)
                    new_cpd = Compound.Compound(new_mol)
                    new_cpd.comments.append(
                        "%d AromaticRingRotate6 %d %d %d %d %d %d" % (i, seq[0], seq[1], seq[2], seq[3], seq[4],
                                                                        seq[5]))
                except:
                    print(sys.exc_info())
                    continue
                if new_inchi in self.inchis: 
                    print("new_inchi in self.inchis")
                    #continue
                else:
                    print("new_inchi NOT in self.inchis")
                self.inchis.append(new_inchi)
                self.cpds.append(new_cpd)
                self.cpds[i].children.append(len(self.cpds) - 1)
                self.cpds[len(self.cpds) - 1].parents.append(i)
        return True

    def AdditionReaction(self, i): # 付加反応 A-B-C <= A=B + C
        for seq in self.cpds[i].FindSeq(3):
            if self.cpds[i].SeqBondOrder(seq, [1, 1]):
                pass
            else:
                continue
            new_cpd = copy.deepcopy(self.cpds[i])
            try:
                new_cpd.AddBond(seq[0] + 1, seq[1] + 1, 2)
                new_cpd.DeleteBond(int(new_cpd.matrix[seq[1]][seq[2]] - 1))
                new_mol = Chem.MolFromMolBlock(new_cpd.GetMolfile())
                new_inchi = Chem.MolToInchi(new_mol)
                new_cpd = Compound.Compound(new_mol)
                new_cpd.comments.append("%d AdditionReaction %d %d %d" % (i, seq[0], seq[1], seq[2]))
            except:
                print(sys.exc_info())
                continue
            if new_inchi in self.inchis: 
                continue
            self.inchis.append(new_inchi)
            self.cpds.append(new_cpd)
            self.cpds[i].children.append(len(self.cpds) - 1)
            self.cpds[len(self.cpds) - 1].parents.append(i)
    
    def EliminationReaction(self, i): # 脱離反応 A=B + C <= A-B-C
        for bond in self.cpds[i].bond_block:
            b1 = bond.split()
            if b1[2] != '2': continue
            for atmidx, atm in enumerate(self.cpds[i].atom_block):
                if atmidx == int(b1[0]) - 1: continue
                if atmidx == int(b1[1]) - 1: continue
                if self.cpds[i].matrix[int(b1[0]) - 1][atmidx] == 0:
                    new_cpd = copy.deepcopy(self.cpds[i])
                    try:
                        new_cpd.AddBond(int(b1[0]), atmidx + 1, 1)
                        new_cpd.AddBond(int(b1[0]), int(b1[1]), 1)
                        new_mol = Chem.MolFromMolBlock(new_cpd.GetMolfile())
                        new_inchi = Chem.MolToInchi(new_mol)
                        new_cpd = Compound.Compound(new_mol)
                        new_cpd.comments.append(
                            "%d EliminationReaction %d %d %d" % (i, int(b1[0]), int(b1[1]), atmidx))
                    except:
                        #print(sys.exc_info())
                        continue
                    if new_inchi in self.inchis: 
                        continue
                    self.inchis.append(new_inchi)
                    self.cpds.append(new_cpd)
                    self.cpds[i].children.append(len(self.cpds) - 1)
                    self.cpds[len(self.cpds) - 1].parents.append(i)

                if self.cpds[i].matrix[int(b1[1]) - 1][i] == 0:
                    new_cpd = copy.deepcopy(self.cpds[i])
                    try:
                        new_cpd.AddBond(int(b1[1]), atmidx + 1, 1)
                        new_cpd.AddBond(int(b1[1]), int(b1[0]), 1)
                        new_mol = Chem.MolFromMolBlock(new_cpd.GetMolfile())
                        new_inchi = Chem.MolToInchi(new_mol)
                        new_cpd = Compound.Compound(new_mol)
                        new_cpd.comments.append(
                            "%d EliminationReaction %d %d %d" % (i, int(b1[1]), int(b1[0]), atmidx))
                    except:
                        #print(sys.exc_info())
                        continue
                    if new_inchi in self.inchis: 
                        continue
                    self.inchis.append(new_inchi)
                    self.cpds.append(new_cpd)
                    self.cpds[i].children.append(len(self.cpds) - 1)
                    self.cpds[len(self.cpds) - 1].parents.append(i)
        return True

    def RearrangementReaction(self, i): # 転位反応 A-B + C <= A + B-C
        for bond1idx, line in enumerate(self.cpds[i].bond_block):
            b1 = line.split()
            if b1[2] != '1': continue
            for atmidx, atm in enumerate(self.cpds[i].atom_block):
                if i == int(b1[0]) - 1: continue
                if i == int(b1[1]) - 1: continue
                if self.cpds[i].matrix[int(b1[0]) - 1][i] == 0:
                    new_cpd = copy.deepcopy(self.cpds[i])
                    try:
                        new_cpd.AddBond(int(b1[0]), atmidx + 1, 1)
                        new_cpd.DeleteBond(bond1idx)
                        new_mol = Chem.MolFromMolBlock(new_cpd.GetMolfile())
                        new_inchi = Chem.MolToInchi(new_mol)
                        new_cpd = Compound.Compound(new_mol)
                        new_cpd.comments.append(
                            "%d RearrangementReaction %d %d %d" % (i, int(b1[0]), int(b1[1]), atmidx))
                    except:
                        #print(sys.exc_info())
                        continue
                    if new_inchi in self.inchis: 
                        continue
                    self.inchis.append(new_inchi)
                    self.cpds.append(new_cpd)
                    self.cpds[i].children.append(len(self.cpds) - 1)
                    self.cpds[len(self.cpds) - 1].parents.append(i)
                if self.cpds[i].matrix[int(b1[1]) - 1][i] == 0:
                    new_cpd = copy.deepcopy(self.cpds[i])
                    try:
                        new_cpd.AddBond(int(b1[1]), atmidx + 1, 1)
                        new_cpd.DeleteBond(bond1idx)
                        new_mol = Chem.MolFromMolBlock(new_cpd.GetMolfile())
                        new_inchi = Chem.MolToInchi(new_mol)
                        new_cpd = Compound.Compound(new_mol)
                        new_cpd.comments.append(
                            "%d RearrangementReaction %d %d %d" % (i, int(b1[1]), int(b1[0]), atmidx))
                    except:
                        #print(sys.exc_info())
                        continue
                    if new_inchi in self.inchis: 
                        continue
                    self.inchis.append(new_inchi)
                    self.cpds.append(new_cpd)
                    self.cpds[i].children.append(len(self.cpds) - 1)
                    self.cpds[len(self.cpds) - 1].parents.append(i)
        return True

    def DiOxidativeDigestion(self, i):
        bonds = []
        for bond1idx, bond in enumerate(self.cpds[i].bond_block):
            b1 = bond.split()
            if b1[2] == '2':
                a1 = self.cpds[i].atom_symbols[int(b1[0]) - 1]
                a2 = self.cpds[i].atom_symbols[int(b1[1]) - 1]
                if a1 == "O" and a2 != "O":
                    bonds.append([bond1idx, int(b1[0]), int(b1[1])])
                elif a1 != "O" and a2 == "O":
                    bonds.append([bond1idx, int(b1[1]), int(b1[0])])
        for ary1 in bonds:
            for ary2 in bonds:
                if ary1[0] >= ary2[0]: continue
                new_cpd = copy.deepcopy(self.cpds[i])
                try:
                    new_cpd.DeleteBond(ary1[0])
                    new_cpd.DeleteBond(ary2[0])
                    new_cpd.AddBond(ary1[1], ary2[1], 2)
                    new_cpd.AddBond(ary1[2], ary2[2], 2)
                    new_mol = Chem.MolFromMolBlock(new_cpd.GetMolfile())
                    new_inchi = Chem.MolToInchi(new_mol)
                    new_cpd = Compound.Compound(new_mol)
                    new_cpd.comments.append(
                        "%d DiOxidativeDigestion %d %d" % (i, ary1[0], ary2[0]))
                except:
                    print(sys.exc_info())
                    continue
                if new_inchi in self.inchis: 
                    continue
                self.inchis.append(new_inchi)
                self.cpds.append(new_cpd)
                self.cpds[i].children.append(len(self.cpds) - 1)
                self.cpds[len(self.cpds) - 1].parents.append(i)
        return True

    def InsertOxygen(self, i): # 酸素挿入 A-O-B <= A-B + O
        seqs = []
        for seq in self.cpds[i].FindUniqSeq(3):
            if self.cpds[i].atom_symbols[seq[1]] != 'O': continue
            new_cpd = copy.deepcopy(self.cpds[i])
            try:
                new_cpd.DeleteBond(new_cpd.matrix[seq[0]][seq[1]] - 1)
                new_cpd.DeleteBond(new_cpd.matrix[seq[1]][seq[2]] - 1)
                new_cpd.AddBond(seq[0] + 1, seq[2] + 1, 1)
                new_mol = Chem.MolFromMolBlock(new_cpd.GetMolfile())
                new_inchi = Chem.MolToInchi(new_mol)
                new_cpd = Compound.Compound(new_mol)
                new_cpd.comments.append(
                    "%d InsertOxygen %d %d %d" % (i, seq[0], seq[1], seq[2]))
            except:
                print(sys.exc_info())
                continue
            if new_inchi in self.inchis: 
                continue
            self.inchis.append(new_inchi)
            self.cpds.append(new_cpd)
            self.cpds[i].children.append(len(self.cpds) - 1)
            self.cpds[len(self.cpds) - 1].parents.append(i)
        return True

    def DielsAlder(self, i): # ディールス・アルダー反応 C1-C=C-C-C-C1 <= C=C-C=C + C=C
        for seq in self.cpds[i].FindUniqRing(6):
            if self.cpds[i].RingBondOrder(seq, [1, 2, 1, 1, 1, 1]):
                new_cpd = copy.deepcopy(self.cpds[i])
                try:
                    new_cpd.AddBond(seq[0] + 1, seq[1] + 1, 2)
                    new_cpd.AddBond(seq[1] + 1, seq[2] + 1, 1)
                    new_cpd.AddBond(seq[2] + 1, seq[3] + 1, 2)
                    new_cpd.AddBond(seq[4] + 1, seq[5] + 1, 2)
                    new_cpd.DeleteBond(int(new_cpd.matrix[seq[3]][seq[4]]) - 1)
                    new_cpd.DeleteBond(int(new_cpd.matrix[seq[0]][seq[5]]) - 1)
                    new_mol = Chem.MolFromMolBlock(new_cpd.GetMolfile())
                    new_inchi = Chem.MolToInchi(new_mol)
                    new_cpd = Compound.Compound(new_mol)
                    new_cpd.comments.append(
                    "%d DielsAlder %d %d %d %d %d %d" % (i, seq[0], seq[1], seq[2], seq[3], seq[4], seq[5]))
                except:
                    print(sys.exc_info())
                    continue
                if new_inchi in self.inchis: 
                    continue
                self.inchis.append(new_inchi)
                self.cpds.append(new_cpd)
                self.cpds[i].children.append(len(self.cpds) - 1)
                self.cpds[len(self.cpds) - 1].parents.append(i)
        return True

    def RetroDielsAlder(self, i): # 逆ディールス・アルダー反応 C=C-C=C + C=C <= C1-C=C-C-C-C1
        #print("RetroDielsAlder")
        for seq in self.cpds[i].FindSeq(4):
            if self.cpds[i].SeqBondOrder(seq, [2, 1, 2]):
                #print(seq)
                pass
            else:
                continue
            for bondidx, bond in enumerate(self.cpds[i].bond_block):
                b4 = bond.split()
                #print(seq, b4)
                if b4[2] != '2': continue
                a1, a2 = int(b4[0]) - 1, int(b4[1]) - 1
                if (a1 in seq) or (a2 in seq): continue
                #print(seq, a1, a2)
                for bb in [[a1, a2], [a2, a1]]:
                    new_cpd = copy.deepcopy(self.cpds[i])
                    try:
                        new_cpd.AddBond(seq[0] + 1, seq[1] + 1, 1)
                        new_cpd.AddBond(seq[1] + 1, seq[2] + 1, 2)
                        new_cpd.AddBond(seq[2] + 1, seq[3] + 1, 1)
                        new_cpd.AddBond(seq[3] + 1, bb[0] + 1, 1)
                        new_cpd.AddBond(bb[0] + 1, bb[1] + 1, 1)
                        new_cpd.AddBond(bb[1] + 1, seq[0] + 1, 1)
                        new_mol = Chem.MolFromMolBlock(new_cpd.GetMolfile())
                        new_inchi = Chem.MolToInchi(new_mol)
                        new_cpd = Compound.Compound(new_mol)
                        new_cpd.comments.append(
                        "%d RetroDielsAlder %d %d %d %d %d %d" % (i, seq[0], seq[1], seq[2], seq[3], bb[0], bb[1]))
                    except:
                        print(sys.exc_info())
                        continue
                    if new_inchi in self.inchis: 
                        continue
                    self.inchis.append(new_inchi)
                    self.cpds.append(new_cpd)
                    self.cpds[i].children.append(len(self.cpds) - 1)
                    self.cpds[len(self.cpds) - 1].parents.append(i)
        return True
  
    def CopeRearrangement(self, i): # コープ転位・クライゼン転位 A=B-C-D-E=F <= C=B-A-F-E=D
        for seq in self.cpds[i].FindSeq(6):
            if self.cpds[i].SeqBondOrder(seq, [2, 1, 1, 1, 2]):
                pass
            else:
                continue
            new_cpd = copy.deepcopy(self.cpds[i])
            try:
                new_cpd.AddBond(seq[0] + 1, seq[1] + 1, 1) # AB
                new_cpd.AddBond(seq[1] + 1, seq[2] + 1, 2) # BC
                new_cpd.DeleteBond(int(new_cpd.matrix[seq[2]][seq[3]]) - 1) # CD
                new_cpd.AddBond(seq[3] + 1, seq[4] + 1, 2) # DE
                new_cpd.AddBond(seq[4] + 1, seq[5] + 1, 1) # EF
                new_cpd.AddBond(seq[0] + 1, seq[5] + 1, 1) # AF
                new_mol = Chem.MolFromMolBlock(new_cpd.GetMolfile())
                new_inchi = Chem.MolToInchi(new_mol)
                new_cpd = Compound.Compound(new_mol)
                new_cpd.comments.append(
                "%d CopeRearrangement %d %d %d %d %d %d" % (i, seq[0], seq[1], seq[2], seq[3], seq[4], seq[5]))
            except:
                print(sys.exc_info())
                continue
            if new_inchi in self.inchis: 
                continue
            self.inchis.append(new_inchi)
            self.cpds.append(new_cpd)
            self.cpds[i].children.append(len(self.cpds) - 1)
            self.cpds[len(self.cpds) - 1].parents.append(i)
        return True

    def Cycloaddition22(self, i): # [2,2]環化付加 A1-B-C-D1 <= A=B + C=D
        for seq in self.cpds[i].FindRing(4):
            if self.cpds[i].RingBondOrder(seq, [1, 1, 1, 1]):
                pass
            else:
                continue
            new_cpd = copy.deepcopy(self.cpds[i])
            try:
                new_cpd.AddBond(seq[0] + 1, seq[1] + 1, 2) # AB
                new_cpd.DeleteBond(int(new_cpd.matrix[seq[1]][seq[2]]) - 1) # BC
                new_cpd.AddBond(seq[2] + 1, seq[3] + 1, 2) # CD
                new_cpd.DeleteBond(int(new_cpd.matrix[seq[3]][seq[0]]) - 1) # DA
                new_mol = Chem.MolFromMolBlock(new_cpd.GetMolfile())
                new_inchi = Chem.MolToInchi(new_mol)
                new_cpd = Compound.Compound(new_mol)
                new_cpd.comments.append(
                "%d Cycloaddition22 %d %d %d %d" % (i, seq[0], seq[1], seq[2], seq[3]))
            except:
                print(sys.exc_info())
                continue
            if new_inchi in self.inchis: 
                continue
            self.inchis.append(new_inchi)
            self.cpds.append(new_cpd)
            self.cpds[i].children.append(len(self.cpds) - 1)
            self.cpds[len(self.cpds) - 1].parents.append(i)
        return True