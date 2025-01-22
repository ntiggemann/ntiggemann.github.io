from kivy.app import App
from kivy.uix.boxlayout import BoxLayout
from kivy.uix.button import Button
from kivy.uix.textinput import TextInput
from kivy.uix.gridlayout import GridLayout
from kivy_garden.matplotlib import FigureCanvasKivyAgg
from kivy.clock import Clock
import matplotlib.pyplot as plt
from kivy.uix.label import Label
import functions_tpa # type: ignore

class PlotTrop(BoxLayout):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.orientation = "vertical"

        # Plot-Anzeige
        self.plot_area = BoxLayout(size_hint_y=2)
        self.plot_area.orientation = "vertical"
        self.add_widget(self.plot_area)

        # Layout für die Grenzen
        bounds_layout = BoxLayout(orientation="horizontal", spacing=10, size_hint=(1, None), height=40)
        self.xmin_input = TextInput(text="-5",hint_text="xmin", multiline=False, size_hint=(1, None), height=40)
        self.xmax_input = TextInput(text="5",hint_text="xmax", multiline=False, size_hint=(1, None), height=40)
        self.ymin_input = TextInput(text="-5",hint_text="ymin", multiline=False, size_hint=(1, None), height=40)
        self.ymax_input = TextInput(text="5",hint_text="ymax", multiline=False, size_hint=(1, None), height=40)

        bounds_layout.add_widget(Label(text="xmin:"))
        bounds_layout.add_widget(self.xmin_input)
        bounds_layout.add_widget(Label(text="xmax:"))
        bounds_layout.add_widget(self.xmax_input)
        bounds_layout.add_widget(Label(text="ymin:"))
        bounds_layout.add_widget(self.ymin_input)
        bounds_layout.add_widget(Label(text="ymax:"))
        bounds_layout.add_widget(self.ymax_input)

        self.add_widget(bounds_layout)

        # Eingabefeld
        self.input_field = TextInput(
            hint_text="Input e.g. -2x2-2y2-2+x+y+xy",
            multiline=False,
            readonly=True,
            size_hint=(1, None),  # Breite flexibel, Höhe fix
            height=40  # Feste Höhe in Pixeln
        )
        self.input_field.bind(focus=Clock.schedule_once(lambda dt: self.toggle_keyboard, 0))
        self.add_widget(self.input_field)

        # Tastatur-Layout
        self.keyboard_layout = GridLayout(cols=3, spacing=3, padding=3, size_hint_y=None)
        self.keyboard_layout.height = 200
        self.keyboard_layout.visible = False  # Anfangszustand: Tastatur unsichtbar

        keys = [
            "7", "8", "9",
            "4", "5", "6",
            "1", "2", "3",
            "+", "0", "-",
            "x", "y", "delete"
        ]

        for key in keys:
            button = Button(
                text=key,
                font_size=24,
                on_press=self.on_key_press,
            )
            self.keyboard_layout.add_widget(button)

        self.add_widget(self.keyboard_layout)

        # Plot-Button
        self.plot_button = Button(
            text="plot curve",
            size_hint_y=None,
            height=50,
            font_size=24,
            on_press=self.plot_polynomial
        )
        self.add_widget(self.plot_button)

    def toggle_keyboard(self, instance, value):
        if value:
            self.keyboard_layout.opacity = 1  # Tastatur anzeigen
        else:
            self.keyboard_layout.opacity = 0  # Tastatur ausblenden

    def on_key_press(self, instance):
        cursor_pos = self.input_field.cursor_index()

        if instance.text == "delete":
            if cursor_pos > 0:
                self.input_field.text = (
                    self.input_field.text[:cursor_pos - 1]
                    + self.input_field.text[cursor_pos:]
                )
                self.input_field.cursor = (cursor_pos - 1, 0)
        else:
            self.input_field.text = (
                self.input_field.text[:cursor_pos]
                + instance.text
                + self.input_field.text[cursor_pos:]
            )
            self.input_field.cursor = (cursor_pos + 1, 0)

    def plot_polynomial(self, instance):
        poly_str = self.input_field.text
        a_x = float(self.xmin_input.text) if self.xmin_input.text else -5
        b_x = float(self.xmax_input.text) if self.xmax_input.text else 5
        a_y = float(self.ymin_input.text) if self.ymin_input.text else -5
        b_y = float(self.ymax_input.text) if self.ymax_input.text else 5
        weight_text_size = 10
        weight_color = "blue"
        color_curve = "blue"
        color_DS = "blue"
        color_NP = "black"
        eps = 0.00000000000001
        try:
            functions_tpa.plot_trop_polynomial(poly_str,a_x,b_x,a_y,b_y,weight_text_size,weight_color,
                        color_curve,color_DS,color_NP,eps)

            # Plot in Kivy anzeigen
            self.plot_area.clear_widgets()
            self.plot_area.add_widget(FigureCanvasKivyAgg(plt.figure(1)))
            self.plot_area.add_widget(FigureCanvasKivyAgg(plt.figure(2)))

        except Exception as e:
            self.input_field.text = "Fehler"

class KeyboardApp(App):
    def build(self):
        return PlotTrop()

if __name__ == "__main__":
    KeyboardApp().run()
