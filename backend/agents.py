import os
from dotenv import load_dotenv
from crewai import Agent, LLM
from backend.tools import search_pubmed

# Load environment variables
load_dotenv()

# Configuration for the Groq LLM
medical_llm = LLM(
    model="groq/llama-3.3-70b-versatile",
    api_key=os.getenv("GROQ_API_KEY"),
    temperature=0.1  # Set low for higher clinical accuracy and relevance
)

# Agent 1: Specialist to refine raw symptoms into clinical search terms
classifier_agent = Agent(
    role="Medical Symptom Triage Specialist",
    goal="Extract and refine medical keywords from the user symptoms: {symptoms}",
    backstory=(
        "You are an expert at medical intake. You convert casual descriptions into "
        "professional clinical MeSH terms to ensure high-quality research matches."
    ),
    llm=medical_llm,
    allow_delegation=False,
    verbose=True
)

# Agent 2: Researcher to find evidence-based papers on PubMed
matcher_agent = Agent(
    role="Clinical Research Librarian",
    goal="Use PubMed to find relevant clinical trials and research papers for {symptoms}.",
    backstory=(
        "You are a master of medical databases. You prioritize human clinical trials "
        "and systematic reviews. You ignore rare case reports unless they are directly "
        "relevant to the user's symptoms to ensure the most helpful results."
    ),
    tools=[search_pubmed],
    llm=medical_llm,
    allow_delegation=False,
    verbose=True
)

# Agent 3: Advisor to compile the final professional report
advisor_agent = Agent(
    role="Medical Research Advisor",
    goal="Synthesize research findings into a structured, professional medical report.",
    backstory=(
        "You translate complex research into clear, actionable information. You ensure "
        "the structure is professional and strictly adheres to safety disclaimers."
    ),
    llm=medical_llm,
    allow_delegation=False,
    verbose=True
)