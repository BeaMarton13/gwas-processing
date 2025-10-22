import pandas as pd
from collections import Counter
import re
import sys 


class ExcelHandler:
    def __init__(self, file_path):
        self.file_path = file_path
        self.df = self.read_excel(file_path)

    
    # def get_disease_column(self, column_name, keyword):
    #     mask = self.df[column_name].notna()
    #     return self.df[mask & self.df[column_name].str.contains(keyword, case=False)]

    def get_disease_column(self, column_name, keyword):
        mask = self.df[column_name].notna()
        pattern = fr"\b{keyword}\b"   # \b ensures whole word match
        return self.df[mask & self.df[column_name].str.contains(pattern, case=False, regex=True)]


    def get_full_ranking(self, column_name):
        # Select the disease/trait column
        diseases = self.df[column_name].dropna().astype(str)

        # Count frequencies
        disease_counts = Counter(diseases)

        # Make ranking table
        ranking = pd.DataFrame(disease_counts.most_common(), columns=["Disease", "Count"])
        return ranking
    

    @staticmethod
    def read_excel(file_path):
        df = pd.read_excel(file_path, header=None)  # load without treating any row as header

        # Make the second row (index 1) the new column names
        df.columns = df.iloc[1]

        # Drop the first two rows (old headers)
        df = df.drop([0, 1]).reset_index(drop=True)
        return df


    @staticmethod
    def get_keyword_ranking(text_series):
        # Define a set of common English stopwords
        stopwords = {
            'the', 'of', 'and', 'a', 'to', 'in', 'for', 'on', 'with', 'by', 'an', 'at', 'from',
            'as', 'is', 'are', 'was', 'were', 'be', 'or', 'that', 'this', 'it', 'not', 'which',
            'but', 'has', 'have', 'had', 'their', 'its', 'other', 'also', 'can', 'such', 'may',
            'than', 'these', 'all', 'more', 'one', 'two', 'three', 'four', 'five', 'six', 'seven',
            'eight', 'nine', 'ten', 'no', 'yes', 'do', 'does', 'did', 'so', 'if', 'into', 'out',
            'about', 'up', 'down', 'over', 'under', 'between', 'among', 'after', 'before', 'during',
            'each', 'any', 'some', 'most', 'many', 'much', 'very', 'just', 'like', 'because', 'etc',
            # stopwords from first-hand results
            "disorders", "first", "field", "level", "data", "elsewhere", "ukb", 
            "occurrence", "levels", "non", "due", "classified", "disease", "diseases",
            "cell", "viral",
            "icd10","other","diseases","disease","infection","infections",
            "intoxications","unspecified","disorder","syndrome","type",
            "first","occurrence","ukb","data","field","and","of","to","the","in",
            "count","percentage","mean","due","not","elsewhere","classified","system",
            "acute","chronic","tissue","levels","reticulocyte"
        }

        # Combine all text into one string
        all_text = " ".join(text_series)

        # Tokenize (split into words, lowercase, remove non-letters)
        tokens = re.findall(r'\b[a-zA-Z]+\b', all_text.lower())

        # Remove stopwords
        tokens = [word for word in tokens if word not in stopwords]

        # Count frequencies
        word_counts = Counter(tokens)

        # Convert to DataFrame for ranking
        keyword_ranking = pd.DataFrame(word_counts.most_common(), columns=['Keyword', 'Count'])
        return keyword_ranking


if __name__ == "__main__":
    file_path = sys.path[0] + "/../../data.xlsx"
    excel_handler = ExcelHandler(file_path)

    # Select the disease/trait column
    disease_col = "Reported trait" 

    # Get full ranking of diseases/traits
    ranking = excel_handler.get_full_ranking(disease_col)
    print(ranking.head(20))

    # Get keyword ranking
    pd.set_option('display.max_rows', None)
    keyword_ranking = excel_handler.get_keyword_ranking(excel_handler.df[disease_col].dropna().astype(str))
    print(keyword_ranking) 

    # Filter for mental health related studies
    pd.set_option('display.max_colwidth', None)
    mental_df = excel_handler.get_disease_column(disease_col, " mental ")
    print(mental_df["Study Accession"])
    print("====================")
    print(mental_df[disease_col])
