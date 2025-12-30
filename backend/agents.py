import os
from dotenv import load_dotenv
from crewai import Agent, LLM
from backend.tools import search_pubmed

load_dotenv()

# ============================================================================
# LLM CONFIGURATION - OPTIMIZED FOR GROQ FREE TIER
# ============================================================================

# THE WORKHORSE: Used for keyword extraction and data searching.
# 8B is fast and has higher rate limits.
fast_llm = LLM(
    model="groq/llama-3.1-8b-instant",
    api_key=os.getenv("GROQ_API_KEY"),
    temperature=0.1,
    max_tokens=1024  # REDUCED from 2048 to prevent overflow
)

# THE BRAIN: Used only for the final synthesis.
# 70B is smart but has very strict token-per-minute limits.
clinical_llm = LLM(
    model="groq/llama-3.3-70b-versatile",
    api_key=os.getenv("GROQ_API_KEY"),
    temperature=0.1,
    max_tokens=2048  # Keep higher for final report
)

# ============================================================================
# AGENT DEFINITIONS - MINIMAL ITERATIONS TO PREVENT TOKEN LOOPS
# ============================================================================

# Agent 1: Extracts clinical terms
classifier_agent = Agent(
    role="Medical Symptom Triage Specialist",
    goal="Extract exactly 2-3 medical keywords from: {symptoms}",
    backstory=(
        "You are an expert at converting patient descriptions into "
        "concise MeSH terms. Keep responses under 20 words."
    ),
    llm=fast_llm,
    max_iter=1,        # CRITICAL: Only 1 iteration to prevent looping
    max_rpm=10,        
    allow_delegation=False,
    verbose=True
)

# Agent 2: Searches and reads PubMed
# CRITICAL FIX: max_iter=1 prevents the agent from retrying on token errors
matcher_agent = Agent(
    role="Clinical Research Librarian",
    goal="Find ONE relevant clinical study and provide a 2-sentence summary.",
    backstory=(
        "You are a master of medical databases. You extract ONLY the title "
        "and conclusion from PubMed studies. Keep all responses under 100 words."
    ),
    tools=[search_pubmed],
    llm=fast_llm, 
    max_iter=1,        # CRITICAL: Only 1 attempt - if tool returns data, use it immediately
    max_rpm=2,         
    allow_delegation=False,
    verbose=True
)

# Agent 3: The Advisor (The "Brain")
# This agent now gets minimal input from Agent 2, so it won't overflow
advisor_agent = Agent(
    role="Medical Research Advisor",
    goal="Create a brief, professional medical research summary in Markdown.",
    backstory=(
        "You are a senior clinical consultant. You take short research summaries "
        "and format them into clear, concise reports with disclaimers. "
        "Keep reports under 400 words."
    ),
    llm=clinical_llm,
    max_iter=1,        # Single iteration for final output
    max_rpm=1,         
    allow_delegation=False,
    verbose=True
)