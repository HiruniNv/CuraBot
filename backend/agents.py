import os
from dotenv import load_dotenv
from crewai import Agent, LLM
from backend.tools import search_pubmed

# Load environment variables from .env
load_dotenv()

# Define the LLM configuration for Groq
# We use Llama-3.3-70b-versatile for high-quality medical reasoning
medical_llm = LLM(
    model="groq/llama-3.3-70b-versatile",
    api_key=os.getenv("GROQ_API_KEY"),
    temperature=0.2 # Lower temperature for higher clinical accuracy
)

# --- Define Agents ---

classifier_agent = Agent(
    role="Medical Symptom Triage Specialist",
    goal="Extract and refine medical keywords from the user symptoms: {symptoms}",
    backstory=(
        "You are an expert at medical intake. You turn casual descriptions into "
        "professional clinical terms (e.g., 'stomach pain' to 'abdominal distress') "
        "to ensure high-quality research matches."
    ),
    llm=medical_llm,
    allow_delegation=False,
    verbose=True
)

matcher_agent = Agent(
    role="Clinical Research Librarian",
    goal="Use PubMed to find the most relevant clinical trials and research papers.",
    backstory=(
        "You are a master of medical databases. You use the search_pubmed tool to "
        "find evidence-based research. You strictly ignore irrelevant geographical results "
        "and focus on pathology and human clinical studies."
    ),
    tools=[search_pubmed],
    llm=medical_llm,
    allow_delegation=False,
    verbose=True
)

advisor_agent = Agent(
    role="Medical Research Advisor",
    goal="Synthesize research findings into a structured, empathetic, and professional medical report.",
    backstory=(
        "You are a senior medical consultant. You translate complex research into "
        "clear, actionable information. You focus on structure, clarity, and "
        "always prioritize patient safety with clear disclaimers."
    ),
    llm=medical_llm,
    allow_delegation=False,
    verbose=True
)