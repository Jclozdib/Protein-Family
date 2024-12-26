import pandas as pd
from wordcloud import WordCloud
import matplotlib.pyplot as plt

file_path = "enrichment_results.tsv"
data = pd.read_csv(file_path, sep="\t")

data["scaled_freq"] = ((1 / data["p_value"]) / (1 / data["p_value"]).max()) * 100

data["scaled_freq"] = data["scaled_freq"].apply(lambda x: max(int(x), 1))

word_freq = {row["name"]: row["scaled_freq"] for _, row in data.iterrows()}

print(word_freq)
wordcloud = WordCloud(
    width=800, height=400, background_color="white"
).generate_from_frequencies(word_freq)

plt.figure(figsize=(10, 5))
plt.imshow(wordcloud, interpolation="bilinear")
plt.axis("off")
plt.show()

output_file = "word_cloud.png"
plt.savefig(output_file)
