import haptk

def read_samples(filename):
    with open(filename) as file:
        lines = [line.rstrip() for line in file]
        return lines

gambian = read_samples("../examples/1kGP_high_coverage_Illumina.gambian.ids")
finnish = read_samples("../examples/1kGP_high_coverage_Illumina.finnish.ids")
han_chinese = read_samples("../examples/1kGP_high_coverage_Illumina.han_chinese.ids")

print(gambian)
print(finnish)
print(han_chinese)

hst = haptk.read_hst("../examples/results/finnish_uhst_left.hst.gz")

hst.circle_tree("fin_my_left_hst.png", to_tag=[gambian, finnish, han_chinese], colors=["red", "blue", "green"])
