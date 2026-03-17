from flask import Flask, jsonify, render_template
import psycopg2
from psycopg2.extras import RealDictCursor

app = Flask(__name__, template_folder="templates")

# ─── DATABASE CONFIG ────────────────────────────────────────────
DB_CONFIG = {
    "host":     "localhost",
    "port":     5432,
    "dbname":   "apo",
    "user":     "song",
    "password": "singSONG!",
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
                    id         SERIAL PRIMARY KEY,
                    hours      INTEGER NOT NULL,
                    pressed_at TIMESTAMPTZ NOT NULL
                );
            """)
        conn.commit()


@app.route("/")
def index():
    return render_template("index.html")


@app.route("/press/<int:hours>", methods=["POST"])
def press(hours):
    if hours not in (0, 1, 2, 3):
        return jsonify({"error": "hours must be 0, 1, 2, or 3"}), 400

    with get_conn() as conn:
        with conn.cursor(cursor_factory=RealDictCursor) as cur:
            cur.execute("""
                INSERT INTO button_presses (hours, pressed_at)
                VALUES (%s, NOW() + %s * INTERVAL '1 hour')
                RETURNING id, hours, pressed_at;
            """, (hours, hours))
            row = cur.fetchone()
        conn.commit()
    return jsonify({
        "id":        row["id"],
        "hours":     row["hours"],
        "pressed_at": row["pressed_at"].isoformat()
    })


@app.route("/history")
def history():
    with get_conn() as conn:
        with conn.cursor(cursor_factory=RealDictCursor) as cur:
            cur.execute("""
                SELECT id, hours, pressed_at
                FROM button_presses
                ORDER BY pressed_at DESC
                LIMIT 10;
            """)
            rows = cur.fetchall()
    return jsonify([
        {"id": r["id"], "hours": r["hours"], "pressed_at": r["pressed_at"].isoformat()}
        for r in rows
    ])


if __name__ == "__main__":
    init_db()
    app.run(host="0.0.0.0", debug=True, port=9000)
