---
outline: deep
---

# Standalone - using Python virtual environment
Cheminformatics Microservice can also be deployed as a standalone application using the Python virtual environment. Here are the step by step instructions.

1. Install Python: Install Python on your machine by following the instructions for your specific operating system.

2. Open a terminal or command prompt.

3. Navigate to the directory where your CM project codebase is located: Use `cd` to navigate to the project directory.

5. Create a virtual environment (optional but recommended): Run the command `python -m venv env` to create a new virtual environment named "env" for your app.

6. Activate the virtual environment (if created): Depending on your operating system, run the appropriate command to activate the virtual environment. For example, on Windows, run `.\env\Scripts\activate`, and on macOS/Linux, run `source env/bin/activate`.

7. Install FastAPI and required dependencies: Run the command `pip install -r requirements.txt` to install FastAPI and the necessary dependencies.

8. Run the FastAPI app: Execute the command `uvicorn main:app --reload` to start the CM app.

9. Wait for the app to start: Uvicorn will start the app and display the server address (usually `http://localhost:8000`) in the terminal or command prompt.

10. Access the FastAPI app: Open a web browser and navigate to the server address displayed in the terminal or command prompt. You should see your FastAPI app running.

That's it!
