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
                    id         INTEGER PRIMARY KEY,
                    pressed_at TIMESTAMPTZ NOT NULL
                );
            """)
        conn.commit()


@app.route("/")
def index():
    return render_template("index.html")


@app.route("/press/<int:hours>", methods=["POST"])
def press(hours):
    if hours not in (1, 2, 3):
        return jsonify({"error": "hours must be 1, 2, or 3"}), 400

    with get_conn() as conn:
        with conn.cursor(cursor_factory=RealDictCursor) as cur:
            cur.execute("""
                INSERT INTO button_presses (id, pressed_at)
                VALUES (1, NOW() + %s * INTERVAL '1 hour')
                ON CONFLICT (id) DO UPDATE
                    SET pressed_at = EXCLUDED.pressed_at
                RETURNING id, pressed_at;
            """, (hours,))
            row = cur.fetchone()
        conn.commit()
    return jsonify({"id": row["id"], "pressed_at": row["pressed_at"].isoformat()})


@app.route("/current")
def current():
    with get_conn() as conn:
        with conn.cursor(cursor_factory=RealDictCursor) as cur:
            cur.execute("SELECT id, pressed_at FROM button_presses WHERE id = 1;")
            row = cur.fetchone()
    if row:
        return jsonify({"pressed_at": row["pressed_at"].isoformat()})
    return jsonify({"pressed_at": None})


if __name__ == "__main__":
    init_db()
    app.run(debug=True, port=9000)
