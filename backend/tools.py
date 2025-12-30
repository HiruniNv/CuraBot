import os
from Bio import Entrez
from crewai.tools import tool 
from dotenv import load_dotenv

load_dotenv()

Entrez.email = os.getenv("NCBI_EMAIL", "your.email@example.com")

@tool("search_pubmed")
def search_pubmed(query: str):
    """
    Searches PubMed for medical literature related to symptoms. 
    The query should be clinical keywords (e.g., 'acute tachycardia human trials').
    Returns a brief summary of the top findings.
    
    CRITICAL: This version is heavily optimized to prevent token overflow errors.
    """
    try:
        # 1. Refine query to pull high-quality evidence only
        # Using [Filter] keeps the results relevant to medical practice
        refined_query = f"{query} AND (clinical trial[Filter] OR review[Filter])"
        
        # 2. Search for IDs - Only get 1 paper to minimize token usage
        handle = Entrez.esearch(db="pubmed", term=refined_query, retmax=1, retmode="json")
        record = Entrez.read(handle)
        handle.close()
        
        id_list = record.get("IdList", [])
        if not id_list:
            return "No high-quality medical literature found for this query. Try broader search terms."

        # 3. Fetch ONLY the title and first 3 sentences of abstract
        # prevents the massive text blocks that crash Groq
        handle = Entrez.efetch(
            db="pubmed", 
            id=id_list[0],  # Only fetch the first result
            rettype="abstract", 
            retmode="xml"  # XML gives more control over parsing
        )
        
        from xml.etree import ElementTree as ET
        tree = ET.parse(handle)
        handle.close()
        root = tree.getroot()
        
        # Extract only what need
        title = root.find(".//ArticleTitle")
        abstract = root.find(".//AbstractText")
        
        title_text = title.text if title is not None else "No title available"
        abstract_text = abstract.text if abstract is not None else "No abstract available"
        
        # --- AGGRESSIVE TRUNCATION FOR TOKEN SAFETY ---
        # Limit abstract to first 500 characters (~125 tokens)
        if len(abstract_text) > 500:
            abstract_text = abstract_text[:500] + "..."
        
        # Construct minimal response
        result = f"**Study Title:** {title_text}\n\n**Key Finding:** {abstract_text}"
        
        # Final safety check: ensure total response is under 800 characters
        if len(result) > 800:
            result = result[:800] + "...[truncated]"
            
        return result

    except Exception as e:
        error_msg = str(e)
        print(f"🔴 PubMed Tool Error: {error_msg}")  # Debug logging
        print(f"🔴 Error Type: {type(e).__name__}")
        
        # Provide helpful error messages
        if "HTTP Error 429" in error_msg or "rate limit" in error_msg.lower():
            return "PubMed rate limit reached. Please wait 1 minute and try again."
        elif "email" in error_msg.lower():
            return "NCBI_EMAIL not configured. Set it in your .env file."
        else:
            return f"PubMed search failed: {error_msg[:100]}"  # Return actual error for debugging