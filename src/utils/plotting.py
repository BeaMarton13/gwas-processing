import matplotlib.pyplot as plt
from adjustText import adjust_text
import numpy as np

def plot_keyword_bubble_chart(keyword_ranking, top_n=30, random_state=42):
    # Select top N keywords
    data = keyword_ranking.head(top_n).copy()
    np.random.seed(random_state)
    # Generate random positions for bubbles
    x = np.random.rand(len(data))
    y = np.random.rand(len(data))
    sizes = np.interp(data['Count'], (data['Count'].min(), data['Count'].max()), (300, 3000))
    labels = [f"{row['Keyword']} ({row['Count']})" for _, row in data.iterrows()]

    fig, ax = plt.subplots(figsize=(12, 8))
    scatter = ax.scatter(x, y, s=sizes, alpha=0.6, c=sizes, cmap='plasma', edgecolors='w', linewidth=2)

    texts = []
    for i, label in enumerate(labels):
        texts.append(ax.text(x[i], y[i], label, ha='center', va='center', fontsize=12, weight='bold'))

    adjust_text(texts, arrowprops=dict(arrowstyle='-', color='gray', lw=0.5))
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title('Keyword Bubble Chart', fontsize=18, weight='bold')
    ax.set_frame_on(False)
    plt.tight_layout()
    plt.show()