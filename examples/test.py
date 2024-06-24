import haptk
# hst = haptk.read_hst("/Users/xosxos/code/haptk/examples/results/finnish_uhst_left.hst.gz")
# hst = haptk.read_hst("/Users/xosxos/code/haptk/examples/results/uhst_left.hst.gz")
hst = haptk.read_hst("/Users/xosxos/code/haptk/examples/results/uhst_right.hst.gz")

def read_samples(filename):
  with open(filename) as file:
      lines = [line.rstrip() for line in file]
      return lines

# The IDs are located in the examples directory
gambian = read_samples("/Users/xosxos/code/haptk/examples/1kGP_high_coverage_Illumina.gambian.ids")
finnish = read_samples("/Users/xosxos/code/haptk/examples/1kGP_high_coverage_Illumina.toscani.ids")
han_chinese = read_samples("/Users/xosxos/code/haptk/examples/1kGP_high_coverage_Illumina.han_chinese.ids")

hst.circle_tree("tmp.png", to_tag=[gambian, finnish, han_chinese], colors=["red", "blue", "green"], w=4560, h=4560)
# hst.circle_tree("tmp.png")
