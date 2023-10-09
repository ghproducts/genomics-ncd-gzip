import pickle
with open('kirnews.pkl', 'rb') as f:
  data = pickle.load(f)

data = pickle.load(open('kirnews.pkl', 'rb'))

print(data['test_data'][0])
print(data['test_label'][0])


#ref = ds['test_labels']
#top_labels = ds['train_labels'][nn]
#n_test = len(ref)

