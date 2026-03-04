# Button Logger

A minimal Flask webapp that logs a UTC timestamp to PostgreSQL each time a button is pressed.

## Setup

### 1. Install dependencies
```bash
pip install -r requirements.txt
```

### 2. Create the Postgres table
Run this once in your database:
```sql
CREATE TABLE IF NOT EXISTS button_presses (
    id         SERIAL PRIMARY KEY,
    pressed_at TIMESTAMPTZ NOT NULL DEFAULT NOW()
);
```

### 3. Set your database credentials
Open `app.py` and fill in the `DB_CONFIG` block near the top:
```python
DB_CONFIG = {
    "host":     "YOUR_HOST",       # e.g. "localhost" or an RDS endpoint
    "port":     5432,
    "dbname":   "YOUR_DATABASE",
    "user":     "YOUR_USER",
    "password": "YOUR_PASSWORD",
}
```

### 4. Run the app
```bash
python app.py
```

Visit http://localhost:5000 — press the button, watch timestamps appear.

## Project layout
```
button-logger/
├── app.py              # Flask backend
├── requirements.txt
└── templates/
    └── index.html      # Frontend UI
```

## API endpoints
| Method | Path       | Description                          |
|--------|------------|--------------------------------------|
| GET    | `/`        | Renders the UI                       |
| POST   | `/press`   | Inserts a timestamp, returns the row |
| GET    | `/history` | Returns last 50 presses as JSON      |
