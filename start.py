#!/usr/bin/env python3
import webbrowser
import threading
import sys
import os
import time

# Ensure we can import app from the current directory
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from app import app

def open_browser():
    """Opens the browser after a short delay to allow the server to start."""
    time.sleep(1.5)
    url = "http://127.0.0.1:5000"
    print(f"Opening browser at {url}...")
    webbrowser.open(url)

if __name__ == "__main__":
    print("Starting Enzyme-Substrate Simulator...")
    
    # Schedule browser launch
    threading.Thread(target=open_browser).start()
    
    # Run the Flask app
    # debug=False to avoid reloader (which might open browser twice or cause issues in a launcher script)
    app.run(port=5000, debug=False)
