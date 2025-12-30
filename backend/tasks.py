from crewai import Task

def create_tasks(user_input, classifier, matcher, advisor):
    """
    Creates a token-optimized task chain for medical research.
    
    KEY OPTIMIZATION: All tasks are deliberately short and direct
    to prevent token overflow errors with Groq free tier.
    """
    
    # Task 1: Extract keywords ONLY - no explanations
    task1 = Task(
        description=(
            f"Input: '{user_input}'\n"
            "Task: List 2-3 MeSH medical terms. Nothing else."
        ),
        expected_output="3 keywords separated by commas.",
        agent=classifier
    )
    
    # Task 2: Search and summarize in ONE sentence
    # CRITICAL: explicitly limit output length to prevent token overflow
    task2 = Task(
        description=(
            "Use the keywords to search PubMed. "
            "Return: Title of 1 study + 1 sentence summary. "
            "Maximum 50 words total."
        ),
        expected_output="Study title and 1-sentence finding.",
        agent=matcher,
        context=[task1]
    )
    
    # Task 3: Final report with strict length limits
    task3 = Task(
        description=(
            "Create a medical report using this template:\n\n"
            "### 🩺 RESEARCH SUMMARY\n"
            "**Input:** [1 sentence]\n\n"
            "### 🔍 FINDING\n"
            "[2 sentences from research]\n\n"
            "### 💡 CONSIDERATIONS\n"
            "[2-3 bullet points]\n\n"
            "### 📋 NEXT STEPS\n"
            "⚠️ This is research info, not diagnosis. See a doctor.\n"
            "[Suggest 1 medical department]\n\n"
            "IMPORTANT: Total report must be under 300 words."
        ),
        expected_output="A complete Markdown report under 300 words.",
        agent=advisor,
        context=[task2]
    )
    
    return [task1, task2, task3]