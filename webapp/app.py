from flask import Flask, jsonify, render_template
import psycopg2
from psycopg2.extras import RealDictCursor
from datetime import datetime

app = Flask(__name__)

# ─── DATABASE CONFIG ────────────────────────────────────────────
DB_CONFIG = {
    "host":     "YOUR_HOST",
    "port":     5432,
    "dbname":   "YOUR_DATABASE",
    "user":     "YOUR_USER",
    "password": "YOUR_PASSWORD",
}
# ────────────────────────────────────────────────────────────────

def get_conn():
    return psycopg2.connect(**DB_CONFIG)


def init_db():
    """Create the button_presses table if it doesn't exist."""
    with get_conn() as conn:
        with conn.cursor() as cur:
            cur.execute("""
                CREATE TABLE IF NOT EXISTS button_presses (
                    id        SERIAL PRIMARY KEY,
                    pressed_at TIMESTAMPTZ NOT NULL DEFAULT NOW()
                );
            """)
        conn.commit()


@app.route("/")
def index():
    return render_template("index.html")


@app.route("/press", methods=["POST"])
def press():
    with get_conn() as conn:
        with conn.cursor(cursor_factory=RealDictCursor) as cur:
            cur.execute(
                "INSERT INTO button_presses (pressed_at) VALUES (NOW()) RETURNING id, pressed_at;"
            )
            row = cur.fetchone()
        conn.commit()
    return jsonify({"id": row["id"], "pressed_at": row["pressed_at"].isoformat()})


@app.route("/history")
def history():
    with get_conn() as conn:
        with conn.cursor(cursor_factory=RealDictCursor) as cur:
            cur.execute(
                "SELECT id, pressed_at FROM button_presses ORDER BY pressed_at DESC LIMIT 50;"
            )
            rows = cur.fetchall()
    return jsonify([
        {"id": r["id"], "pressed_at": r["pressed_at"].isoformat()}
        for r in rows
    ])


if __name__ == "__main__":
    init_db()
    app.run(debug=True, port=5000)
