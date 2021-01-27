#入力した組成式を分解して、画像形式へ。（族,周期）形式で、族は32個で、周期は7
from pymatgen.core import Composition
from pymatgen.core import periodic_table
import json
import matplotlib.pyplot as plt
import numpy as np
import re
from matplotlib.colors import LogNorm

#元素を番号順にソートしlistで出力
# How to use: print(get_elements_list()) >> ['Non','H','He',........'Lr','Non'...] total length is 103(Rf to Og is removed)
def get_elements_list():    
    with open('./periodic_table.json') as f:
        periodic_table_json = json.load(f)
    #各元素の番号を抜き出し、その要素にイオンをぶち込み
    elements_list = ["Non"] * 118
    for el in periodic_table_json.keys():
        el_num = periodic_table_json[el]["Atomic no"]
        elements_list[el_num] = el
    return elements_list

#入力元素の元素番号を取得
# How to use: print(get_element_num("Pb")) >> 82
def get_element_num(element):
    with open('./periodic_table.json') as f:
        periodic_table_json = json.load(f)
    return periodic_table_json[element]["Atomic no"]

#入力組成を分解するクラス
# How to use
# decomposed_comp = comp_decompose("PbMnO3")
# print(decomposed_comp.get_dict()) >> {'Pb': 0.2, 'Mn': 0.2, 'O': 0.6}
# print(decomposed_comp.get_ratio_list()) >> [0,0,0,...0.6,0,0...0.2,....0,0.2.....0] 周期表に沿った構成比を出力
class comp_decompose:
    ratio_list = [0] * 118 #空の元素比率リスト
    def __init__(self, composition):
        self.comp = Composition(composition)
        
    def get_dict(self):
        return self.comp.fractional_composition.as_dict()
    def get_ratio_list(self):
        for el,ratio in self.comp.fractional_composition.as_dict().items():
            self.ratio_list[get_element_num(el)] = ratio
        return self.ratio_list
    def get_normalized_magnification(self):
        comp = self.comp
        parent_ratio = [ratio[1] for ratio in comp.as_dict().items()]
        converted_ratio = [fratio[1] for fratio in comp.fractional_composition.as_dict().items()]
        normalized_magnification = parent_ratio[0]/converted_ratio[0]
        return normalized_magnification

#拡張周期表の元素の族と周期を計算
#拡張周期表とは、ランタノイドなども列として加えたもので、画像としてインプットするための変換である。Lnは15種あるため、例えばHeは1周期の17+15=32属となる。
# How to use
# element_by_PT = element_format_to_PT("Pb") 
# print(element_by_PT.get_row()) >> 6 (６周期)
# print(element_by_PT.get_group()) >> 28(拡張周期表における28属)
# print(element_by_PT.get_converted_number()) >> 188 (拡張周期表で188番目の元素)
class element_format_to_PT:
    def __init__(self, element):
        self.element = periodic_table.Element(element)
        #元素の族と周期から拡張周期表の族と周期に変換
        group = self.element.group
        row = self.element.row
        if group > 3:
            if row < 8:
                self.converted_row = row
                self.converted_group = group + 14
            else:
                self.converted_row = row -2
                self.converted_group = group
        else:
            if row < 8:
                self.converted_row = row
                self.converted_group = group
            else:
                self.converted_row = row -2
                self.converted_group = group
    def get_row(self):#周期を取得
        return self.converted_row
    def get_group(self):#族を取得
        return self.converted_group
    def get_converted_number(self):#拡張周期表の何番目かを取得
        return (self.converted_row-1)*32 + self.converted_group -1 #Hを左上にするための-1

#組成から拡張周期表のヒートマップ画像等を生成
# How to use
# comp_PT_imger = generate_extended_PT_from_composition("PbMnO3")
# print(comp_PT_imger.get_1d_array()) >> [0,0,0,......0.6,......0.2,0,.....0.2,0,...] length: 224
# print(comp_PT_imger.get_2d_array()) >> [[0,0,0,......],[0,0,...0.6,0,...],[],[],[],[],[]] 各周期ごとにlist化(7x32)
# print(comp_PT_imger.draw_PT_img()) >> 7x32のヒートマップ画像として出力
class generate_extended_PT_from_composition:
    def __init__(self, composition):
        self.comp = composition
        self.has_Os = False
    def get_normalized_magnification(self):
        return comp_decompose(self.comp).get_normalized_magnification()
    def get_1d_array(self):
        PT_extended_ratio = np.zeros(32*7)
        for el,ratio in comp_decompose(self.comp).get_dict().items(): #組成から各元素の名前と比を抽出（e.g. Ba,0.2)
            extended_number = element_format_to_PT(el).get_converted_number() #元素の配列番号を取得
            PT_extended_ratio[extended_number] = ratio #拡張周期表の元素位置に比を挿入
        return PT_extended_ratio
    def get_2d_array(self):
        return self.get_1d_array().reshape((7,32))
    def draw_PT_img(self):
        plt.imshow(self.get_2d_array())
        #,norm=LogNorm()
        #,cmap="gray"
        plt.show()
    def removeOxygen(self):
        #pattern = '(O([0-9\.]*))'
        pattern = '(O[0-9\.]+((\-|\+)[a-zA-Z])?|O((\-|\+)[a-zA-Z])|O[a-rt-zXYZ]?)'
        #pattern = '(O[^[sA-Z]]?[0-9\.]*(\-|\+)?[a-zA-Z]?)'
        #pattern = '(O[0-9\.]*[^s](\-|\+)*[a-zA-Z]*)'
        removed = re.sub(pattern,"",self.comp)
        #print("---" + removed + "---")
        #酸素を切り取った組成をクラス変数へセット
        self.comp = removed
        return self
    def hasOs(self):
        pattern = 'Os'
        result = re.findall(pattern,self.comp)
        if result == []:
          return False
        else:
          return True

