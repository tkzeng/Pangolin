from pkg_resources import resource_filename
from pangolin.model import *

# Change this to the desired models. The model that each number corresponds to is listed below.
model_nums = [0]
# 0 = Heart, P(splice)
# 1 = Heart, usage
# 2 = Liver, P(splice)
# 3 = Liver, usage
# 4 = Brain, P(splice)
# 5 = Brain, usage
# 6 = Testis, P(splice)
# 7 = Testis, usage

# Change this to the desired sequences and strand for each sequence. If the sequence is N bases long, Pangolin will
# return scores for the middle N-10000 bases (so if you are interested in the score for a single site, the input should
# be: 5000 bases before the site, base at the site, 5000 bases after the site). Sequences < 10001 bases can be padded with 'N'.
seqs = [10001*'A']
strands = ['-']

# Load models
models = []
for i in model_nums:
    for j in range(1, 6):
        model = Pangolin(L, W, AR)
        if torch.cuda.is_available():
            model.cuda()
            weights = torch.load(resource_filename("pangolin","models/final.%s.%s.3" % (j, i)))
        else:
            weights = torch.load(resource_filename("pangolin","models/final.%s.%s.3" % (j, i)),
                                 map_location=torch.device('cpu'))
        model.load_state_dict(weights)
        model.eval()
        models.append(model)

# Get scores

IN_MAP = np.asarray([[0, 0, 0, 0],
                     [1, 0, 0, 0],
                     [0, 1, 0, 0],
                     [0, 0, 1, 0],
                     [0, 0, 0, 1]])
INDEX_MAP = {0:1, 1:2, 2:4, 3:5, 4:7, 5:8, 6:10, 7:11}

def one_hot_encode(seq, strand):
    seq = seq.upper().replace('A', '1').replace('C', '2')
    seq = seq.replace('G', '3').replace('T', '4').replace('N', '0')
    if strand == '+':
        seq = np.asarray(list(map(int, list(seq))))
    elif strand == '-':
        seq = np.asarray(list(map(int, list(seq[::-1]))))
        seq = (5 - seq) % 5  # Reverse complement
    return IN_MAP[seq.astype('int8')]

for i, seq in enumerate(seqs):
    seq = one_hot_encode(seq, strands[i]).T
    seq = torch.from_numpy(np.expand_dims(seq, axis=0)).float()

    if torch.cuda.is_available():
        seq = seq.to(torch.device("cuda"))

    for j, model_num in enumerate(model_nums):
        score = []
        # Average across 5 models
        for model in models[5*j:5*j+5]:
            with torch.no_grad():
                score.append(model(seq)[0][INDEX_MAP[model_num],:].cpu().numpy())
        print(np.mean(score, axis=0))
