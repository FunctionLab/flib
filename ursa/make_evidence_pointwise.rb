require 'pp'

def calc_roc(posset, negset)
  total = posset.length * negset.length
  cnt_swap = 0
  cnt_neg = 0

  pos = posset.map{|val| [val,1]}
  neg = negset.map{|val| [val,-1]}

  rank = (pos + neg).sort_by{|tok| - tok[0]}
  rank.each do |r|
    if r[1] == 1
      cnt_swap += cnt_neg
    elsif r[1] == -1
      cnt_neg += 1
    end
  end

  return 1 - (cnt_swap.to_f / total)
end

#pointwise auc
def pointwise(bins, pos_vals, neg_vals)

  pos_hist = {}
  neg_hist = {}
  bins.each do |bin|
    neg_hist[bin] = pos_vals.count{|val| val > bin} + 1 #laplace smoothing
    pos_hist[bin] = neg_vals.count{|val| val < bin} + 1 #laplace smoothing
  end
#pp pos_hist
#pp pos_hist
#exit

  # normalize
  # array.inject(:+)
  neg_norm_fctr = neg_hist.values.inject(:+)
  pos_norm_fctr = pos_hist.values.inject(:+)
  bins.each do |bin|
    neg_hist[bin] = (neg_hist[bin].to_f / neg_norm_fctr).round(8)
    pos_hist[bin] = (pos_hist[bin].to_f / pos_norm_fctr).round(8)
  end

  pos_hist = check_distribution(pos_hist, bins)
  neg_hist = check_distribution(neg_hist, bins)

#pp pos_hist
#gets
#exit
 return pos_hist, neg_hist
end

def histogram(bins, pos_vals, neg_vals)
  pos_hist = distribute(pos_vals, bins)
  neg_hist = distribute(neg_vals, bins)
  return pos_hist, neg_hist
end

def distribute(vals, bins)
  hist = {}
  norm_fctr = bins.length
  bin_size = (bins[1] - bins[0]).round(4) # temporary hack

  ## initialize ; laplace smoothing
  bins.each do |bin|
    hist[bin] ||= 1 
  end

  vals.each do |val|
    #corner case
    max_bin = bins[-1]
    if val > max_bin
      hist[max_bin] += 1
      norm_fctr += 1
      next
    end
    min_bin = bins[0]
    if val <= min_bin
      hist[min_bin] += 1
      norm_fctr += 1
      next
    end

    bins.each do |bin|
      if val > (bin - bin_size).round(4) and val <= bin #original
        hist[bin] += 1
        norm_fctr += 1
        break
      end
    end

  end

  # normalize
  bins.each do |bin|
    #hist[bin] = (hist[bin].to_f / norm_fctr).round(4)
    hist[bin] = (hist[bin].to_f / norm_fctr).round(7)
  end

  #check distribution
  arry = check_distribution(hist, bins)
  arry
end

def check_distribution(hist, bins)
  #check for no zero probabilities
  arry = []
  sum = 0.0
#  while hist.values.inject(:+) > 1.0
#    #pp hist
#    index = bins.reverse.map{|bin| hist[bin]}.each_with_index.max[1]
 #   hist[bins[index]] -= 0.000001
 #   hist[bins[index]] -= 0.0000005
 # end
#pp hist
  bins.each do |bin|
    raise "a bin has zero probability is empty: #{hist}" if hist[bin] == 0.0
    raise "a bin has negative probability: #{bins}" if hist[bin] < 0.0
    raise "a bin has >= 1 probability: #{bins}" if hist[bin] >= 1.0
    #pp hist[bin]
    sum += hist[bin]
    arry << hist[bin]
  end
  # check if sums to 1
#pp sum
  if sum != 1.0
    arry[arry.each_with_index.max[1]] += (1.0 - sum).abs ## some assumptions here..
  end
#pp arry.inject(:+)
  arry
end

def compute_distribution(cv_fn)

  ## parse cv predictions
  lines = File.readlines(cv_fn)
  pos_vals = []
  neg_vals = []
  lines.each do |line|
    toks = line.split("\t")
    if toks[1].to_i == 1
      pos_vals << toks[2].to_f 
    elsif toks[1].to_i == -1
      neg_vals << toks[2].to_f
    end
  end
  pos_vals.sort!
  neg_vals.sort!
 
  num_pos = pos_vals.length
  raise 'number positives < 2' if num_pos < 2
  num_neg = neg_vals.length
  pos_min = pos_vals.min
  pos_max = pos_vals.max
  neg_min = neg_vals.min
  
  ## set bins
#=begin
  bin_size_fctr = Math.sqrt(num_pos).ceil
  bin_size = ((pos_max - pos_min) / bin_size_fctr).round(4)
  
  bins = []
  #bin = neg_vals[num_neg * 1/2] # first bin is the median negative value?
  bin = pos_min.round(4) - 2 * bin_size # first bin is the median negative value?
  bin = neg_min.round(4) - 2 * bin_size # first bin is the median negative value?
  while bin <= pos_max  
    bin = (bin + bin_size).round(4)
    bins << bin
  end
  bins.sort!
#=end

=begin
  vals = (pos_vals+neg_vals).sort
  bins = []
  vals.each_with_index do |val,index|
    bins << val if index % (vals.length/100) == 0
  end
=end
  #bins = (pos_vals+neg_vals).sort.each_with_index.select{|bin, index| index % 100 == 0}[0]
# note vals.length the same across terms
  #  pp bins
#  pp bins.length
#  gets
  ## set probability distribution
  #pos_hist = {} # distribution conditioned that the sample is positive
  #neg_hist = {} # distribution conditioned that the sample is negative

  pos_vals.sort!
  neg_vals.sort!
  #pos_hist, neg_hist = histogram(bins, pos_vals, neg_vals)
  pos_hist, neg_hist = pointwise(bins, pos_vals, neg_vals)

  return [pos_hist, neg_hist, bins]
end


run = 0
(0..0).each do |run|
#out_fn = "holdout/m#{run}-noroot.txt"
out_fn = "holdout/pp#{run}.txt"
out_fn = "all/ppall.txt"

next if File.exist? out_fn

cv_dir = "../../cv/holdout/onto/#{run}/"
cv_dir = "../../cv/all/onto/"

cv_fns = Dir.glob("#{cv_dir}/*.output")
cv_fns.sort!

uid2dstr = {} #mesh uid to histogram distribution
infile = File.open(out_fn,'w')
cv_fns.each do |cv_fn|
  #uid = cv_fn[/D\d+/]
  uid = cv_fn[/BTO\d+/]
#  next if uid.eql? "D000000"
 pp uid
  #puts uid  
  h = compute_distribution(cv_fn)
  #pp h[0]
  #uid2dstr[uid] = h
  infile.puts "#{uid}\t#{h[2].join("|")}\t#{h[0].join("|")}\t#{h[1].join("|")}"
end
infile.close
break
end
