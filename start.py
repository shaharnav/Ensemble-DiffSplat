#!/usr/bin/env python3
import webbrowser
import threading
import sys
import os
import time

# Ensure we can import app from the current directory
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from app import app

import signal

def open_browser():
    """Opens the browser after a short delay to allow the server to start."""
    time.sleep(1.5)
    url = "http://127.0.0.1:5001"
    print(f"Opening browser at {url}...")
    webbrowser.open(url)

def listen_for_close():
    """Listens for terminal input and kills the server gracefully on 'close'."""
    while True:
        try:
            cmd = input("Type 'close' to stop the server and free port 5001: \n")
            if cmd.strip().lower() == "close":
                print("\nShutting down server...")
                os.kill(os.getpid(), signal.SIGINT)
                break
        except (EOFError, KeyboardInterrupt):
            break

if __name__ == "__main__":
    print("Starting Enzyme-Substrate Simulator...")
    
    # Schedule browser launch
    threading.Thread(target=open_browser, daemon=True).start()
    
    # Listen for 'close' command
    threading.Thread(target=listen_for_close, daemon=True).start()
    
    # Run the Flask app
    # debug=False to avoid reloader (which might open browser twice or cause issues in a launcher script)
    app.run(port=5001, debug=False)
