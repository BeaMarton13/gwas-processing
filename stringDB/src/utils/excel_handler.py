"""
Query and analyze the GWAS Catalog spreadsheet (data.xlsx).

Provides utilities to filter studies by trait keywords, rank traits by frequency,
and extract corpus statistics from the "Reported trait" column.
"""
import pandas as pd
from collections import Counter
import re
import sys


class ExcelHandler:
    """Wrapper for querying the GWAS Catalog Excel spreadsheet.

    The GWAS Catalog spreadsheet has metadata in the first row and column headers
    in the second row. This class handles the header offset and provides methods
    to filter studies by keyword and analyze trait text.

    Attributes:
        file_path: Path to the data.xlsx file.
        df: pandas DataFrame of the parsed spreadsheet (header = row 2).
    """
    def __init__(self, file_path):
        """Initialize the ExcelHandler by reading and parsing the Excel file.

        Args:
            file_path: Path to the GWAS Catalog .xlsx file.
        """
        self.file_path = file_path
        self.df = self.read_excel(file_path)

    
    def get_disease_column(self, column_name, keyword, keyword_w_spaces=None, whole_word_match=True):
        """Filter studies where a column contains a keyword (case-insensitive).

        Args:
            column_name: Column to search (e.g., "Reported trait").
            keyword: Keyword for filtering (used if keyword_w_spaces is None).
            keyword_w_spaces: Keyword with spaces (for search). Defaults to keyword.
            whole_word_match: If True, uses regex word boundaries (\\b); otherwise substring match.

        Returns:
            pd.DataFrame: Filtered subset of self.df where the column matches the keyword.
        """
        if keyword_w_spaces is None:
            keyword_w_spaces = keyword
        mask = self.df[column_name].notna()
        if whole_word_match:
            pattern = fr"\b{keyword_w_spaces}\b"   # \b ensures whole word match
        else:
            pattern = keyword_w_spaces
        return self.df[mask & self.df[column_name].str.contains(pattern, case=False, regex=True)]


    def get_full_ranking(self, column_name):
        """Count the frequency of each unique value in a column and return a ranked table.

        Args:
            column_name: Column to analyze (e.g., "Reported trait").

        Returns:
            pd.DataFrame: Ranked table with columns ['Disease', 'Count'], sorted by descending count.
        """
        # Select the disease/trait column
        diseases = self.df[column_name].dropna().astype(str)

        # Count frequencies
        disease_counts = Counter(diseases)

        # Make ranking table
        ranking = pd.DataFrame(disease_counts.most_common(), columns=["Disease", "Count"])
        return ranking
    

    @staticmethod
    def read_excel(file_path):
        """Read the GWAS Catalog Excel file, skipping the metadata row.

        The first row of the file contains metadata; the second row contains
        the actual column headers. This method drops the first two rows and
        uses row 2 (index 1) as the header.

        Args:
            file_path: Path to the .xlsx file.

        Returns:
            pd.DataFrame: Parsed spreadsheet with proper column headers.
        """
        df = pd.read_excel(file_path, header=None)  # load without treating any row as header

        # Make the second row (index 1) the new column names
        df.columns = df.iloc[1]

        # Drop the first two rows (old headers)
        df = df.drop([0, 1]).reset_index(drop=True)
        return df


    @staticmethod
    def get_keyword_ranking(text_series):
        """Tokenize text and rank non-stopword keywords by frequency.

        Combines all text in the series, tokenizes (keeping only alphabetic words),
        removes a large set of English and domain-specific stopwords, and counts
        remaining keyword frequencies.

        Args:
            text_series: pandas Series of strings (e.g., "Reported trait" column).

        Returns:
            pd.DataFrame: Ranked table with columns ['Keyword', 'Count'], sorted by descending count.
        """
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
